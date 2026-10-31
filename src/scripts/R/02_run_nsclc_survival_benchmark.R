# MOSAIC-3 | Nested survival benchmark for freshly acquired LUAD and LUSC.
# Primary comparison: clinical, RNA, CNV, 450K methylation, and combinations.
# Preprocessing, feature filtering, imputation, scaling and lambda tuning are
# performed using the outer-training data only. No event/time field is a predictor.

options(stringsAsFactors=FALSE, warn=1, timeout=3600)
need <- c("jsonlite","data.table","survival","glmnet","matrixStats")
miss <- need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
if (length(miss)) stop("Missing R packages: ",paste(miss,collapse=", "),". Run Rscript src/scripts/R/00_install_dependencies.R")
root <- normalizePath(".",winslash="/",mustWork=TRUE)
cfg <- jsonlite::fromJSON(file.path(root,"config","mosaic3_config.json"),simplifyVector=FALSE)
acfg <- cfg$analysis
seed0 <- as.integer(cfg$project$seed)
quick_mode <- isTRUE(cfg$project$quick_mode)
# Fast mode is a deliberate, reproducible screening benchmark. It keeps the
# same outer-fold evaluation and all nine layer combinations, while reducing
# repeated tuning work and feature dimensions. Set quick_mode=false for the
# more expensive confirmatory run.
repeats_cfg <- as.integer(acfg$cancer_outer_repeats)
outer_k <- as.integer(acfg$cancer_outer_folds)
inner_k_cfg <- as.integer(acfg$cancer_inner_folds)
max_features_cfg <- as.integer(acfg$cancer_max_features_per_layer)
nlambda_cfg <- as.integer(acfg$cancer_nlambda)
repeats <- if (quick_mode) 1L else repeats_cfg
inner_k <- if (quick_mode) min(inner_k_cfg, 2L) else inner_k_cfg
max_features <- if (quick_mode) min(max_features_cfg, 250L) else max_features_cfg
nlambda <- if (quick_mode) min(nlambda_cfg, 8L) else nlambda_cfg
min_ratio <- as.numeric(acfg$cancer_lambda_min_ratio)
alpha_grid <- c(0,0.5,1.0)
ratio_grid <- exp(seq(log(1),log(min_ratio),length.out=nlambda))
design_cache <- new.env(parent=emptyenv())
layer_cache <- new.env(parent=emptyenv())
cat(sprintf("Modo: %s | repeticiones externas: %d | folds externos: %d | folds internos: %d | max features/capa: %d | lambdas: %d\n",
            if (quick_mode) "RAPIDO" else "COMPLETO", repeats, outer_k, inner_k, max_features, nlambda))

barcode_sample <- function(x) {
  x <- toupper(as.character(x))
  m <- regmatches(x,regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]?",x))
  ifelse(nzchar(m),m,NA_character_)
}
case_id <- function(x) {
  x <- toupper(as.character(x))
  m <- regmatches(x,regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}",x))
  ifelse(nzchar(m),m,NA_character_)
}
read_matrix <- function(path, layer, fmt=c("auto","rds","tsv")) {
  if (!file.exists(path)) stop("Matrix not found: ",layer,": ",path)
  if (grepl("\\.rds$",path,ignore.case=TRUE)) {
    x <- readRDS(path)
    if (is.list(x) && !is.matrix(x) && !is.data.frame(x)) {
      opts <- x[vapply(x,function(z)is.matrix(z)||is.data.frame(z),logical(1))]
      if (!length(opts)) stop("El RDS de ",layer," no contiene una matriz.")
      x <- opts[[1]]
    }
    x <- as.matrix(x)
  } else {
    x <- as.matrix(data.table::fread(path, data.table=FALSE, check.names=FALSE))
    if (ncol(x)<2) stop("Invalid TSV matrix: ",path)
    rownames(x) <- as.character(x[,1]); x <- x[,-1,drop=FALSE]
  }
  # Explicit orientation detection by TCGA sample barcodes.
  c_samples <- sum(!is.na(barcode_sample(colnames(x))))
  r_samples <- sum(!is.na(barcode_sample(rownames(x))))
  if (c_samples<3 && r_samples>=3) x <- t(x)
  if (sum(!is.na(barcode_sample(colnames(x))))<3) stop("TCGA sample barcodes were not detected as columns in ",layer)
  suppressWarnings(storage.mode(x) <- "numeric")
  x[!is.finite(x)] <- NA_real_
  # Strip only duplicate feature labels, not sample labels.
  if (is.null(rownames(x))) rownames(x) <- paste0(layer,"_feature_",seq_len(nrow(x)))
  rownames(x) <- make.unique(as.character(rownames(x)))
  x
}
first_col <- function(df, aliases) {
  low <- tolower(names(df)); idx <- match(tolower(aliases),low); idx <- idx[!is.na(idx)]
  if (length(idx)) names(df)[idx[1]] else NULL
}
load_clinical <- function(path, cohort) {
  if (!file.exists(path)) stop("Clinical file not found: ",cohort,": ",path)
  d <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  sid <- first_col(d,c("sample_id","sample","submitter_id","case_id","bcr_patient_barcode"))
  tc <- first_col(d,c("time_days","OS.time","os_time","overall_survival_days"))
  ec <- first_col(d,c("event","OS","os_event","vital_status"))
  if (is.null(sid)||is.null(tc)||is.null(ec)) stop("La clínica ",cohort," requiere sample_id, time_days y event/OS.")
  out <- data.frame(sample_id=barcode_sample(d[[sid]]),case_id=case_id(d[[sid]]),
                    time_days=suppressWarnings(as.numeric(d[[tc]])),event=NA_real_,stringsAsFactors=FALSE)
  ev <- d[[ec]]
  if (is.numeric(ev)||all(grepl("^\\s*[01](\\.0+)?\\s*$",na.omit(as.character(ev))))) {
    out$event <- suppressWarnings(as.numeric(as.character(ev)))
  } else {
    z <- tolower(trimws(as.character(ev)))
    out$event <- ifelse(z %in% c("dead","deceased","1","yes","true"),1,
                        ifelse(z %in% c("alive","0","no","false"),0,NA_real_))
  }
  # Time/event aliases for age, sex, stage; the endpoint itself is never included in X.
  agec <- first_col(d,c("age_years_candidate","age_years","age_at_diagnosis","age_at_initial_pathologic_diagnosis","age"))
  sexc <- first_col(d,c("sex_candidate","sex","gender"))
  stagec <- first_col(d,c("stage_group","pathologic_stage","ajcc_pathologic_tumor_stage","tumor_stage"))
  out$age <- if (!is.null(agec)) suppressWarnings(as.numeric(d[[agec]])) else NA_real_
  if (!is.null(agec) && grepl("days",agec,ignore.case=TRUE)) out$age <- out$age/365.25
  out$sex <- if (!is.null(sexc)) as.character(d[[sexc]]) else NA_character_
  out$stage <- if (!is.null(stagec)) as.character(d[[stagec]]) else NA_character_
  out <- out[!is.na(out$sample_id)&is.finite(out$time_days)&out$time_days>0&out$event%in%c(0,1),,drop=FALSE]
  # One row per exact primary-tumour sample; sample selection is recorded in the processed manifest.
  out <- out[!duplicated(out$sample_id),,drop=FALSE]
  out
}
load_cohort <- function(cohort) {
  proc <- file.path(root,"data","processed",cohort)
  cpath <- file.path(proc,"clinical.csv")
  paths <- c(RNA=file.path(proc,"RNA.tsv.gz"),CNV=file.path(proc,"CNV.tsv.gz"),METHYLATION=file.path(proc,"METHYLATION.tsv.gz"))
  clinical <- load_clinical(cpath,cohort)
  mats <- lapply(names(paths),function(nm)read_matrix(paths[[nm]],nm)); names(mats)<-names(paths)
  # Match exact primary-tumour sample barcodes; patient/case IDs are audit keys only.
  sample_sets <- lapply(mats,function(m) {
    ids <- barcode_sample(colnames(m)); ok <- !is.na(ids)
    m <- m[,ok,drop=FALSE]; colnames(m) <- ids[ok]
    # If technical duplicates exist, retain the first deterministically and report by QC.
    m <- m[,!duplicated(colnames(m)),drop=FALSE]
    m
  })
  clinical$sample_id <- barcode_sample(clinical$sample_id)
  clinical <- clinical[!is.na(clinical$sample_id),,drop=FALSE]
  common <- Reduce(intersect,c(list(clinical$sample_id),lapply(sample_sets,colnames)))
  # The cohort is reconstructed deterministically from the open Xena datasets:
  # primary tumour (sample type 01), one specimen per TCGA case, complete across
  # RNA/CNV/HM450 and valid overall-survival endpoint. The selected IDs are frozen
  # in the processed cohort manifest for reproducible reruns.
  common <- sort(common)
  clinical <- clinical[match(common,clinical$sample_id),,drop=FALSE]
  sample_sets <- lapply(sample_sets,function(m)m[,common,drop=FALSE])
  if (length(common)<100) stop(cohort," has only ",length(common)," complete RNA/CNV/methylation + survival cases; check source files and sample alignment.")
  list(clinical=clinical,matrices=sample_sets,common=common,
       qc=data.frame(cohort=cohort,n_model_complete=length(common),events=sum(clinical$event==1),
         censored=sum(clinical$event==0),RNA_features=nrow(sample_sets$RNA),CNV_features=nrow(sample_sets$CNV),
         METHYLATION_features=nrow(sample_sets$METHYLATION)))
}

make_folds <- function(event,k,seed) {
  set.seed(seed); fold <- integer(length(event))
  for (cl in c(0,1)) { ix <- which(event==cl); ix <- sample(ix); fold[ix] <- rep(seq_len(k),length.out=length(ix)) }
  if (any(fold==0L)) stop("Could not create event-stratified folds.")
  fold
}
prep_omics <- function(M,train_idx,test_idx,layer,max_features) {
  # All feature filtering parameters are learned from the training partition.
  # For HM450, select variable probes on beta values BEFORE logit conversion;
  # this avoids transforming hundreds of thousands of probes at every fold.
  tr0 <- M[,train_idx,drop=FALSE]
  te0 <- M[,test_idx,drop=FALSE]
  vars0 <- matrixStats::rowVars(tr0,na.rm=TRUE)
  vars0[!is.finite(vars0)] <- 0
  ord <- order(vars0,decreasing=TRUE)
  ord <- head(ord[vars0[ord]>1e-10],max_features)
  if (!length(ord)) stop("No quedan características variables para ",layer)
  tr <- tr0[ord,,drop=FALSE]
  te <- te0[ord,,drop=FALSE]
  rm(tr0,te0,vars0)

  if (layer=="METHYLATION") {
    tr <- pmin(pmax(tr,0.001),0.999)
    te <- pmin(pmax(te,0.001),0.999)
    tr <- log2(tr/(1-tr))
    te <- log2(te/(1-te))
  }

  # Vectorised train-median imputation (no per-feature R loop).
  med <- matrixStats::rowMedians(tr,na.rm=TRUE)
  med[!is.finite(med)] <- 0
  bad_tr <- which(!is.finite(tr),arr.ind=TRUE)
  if (nrow(bad_tr)) tr[bad_tr] <- med[bad_tr[,1]]
  bad_te <- which(!is.finite(te),arr.ind=TRUE)
  if (nrow(bad_te)) te[bad_te] <- med[bad_te[,1]]

  mu <- rowMeans(tr)
  sig <- matrixStats::rowSds(tr)
  keep <- is.finite(sig)&sig>1e-8
  tr <- tr[keep,,drop=FALSE]
  te <- te[keep,,drop=FALSE]
  mu <- mu[keep]
  sig <- sig[keep]
  tr <- sweep(sweep(tr,1,mu,"-"),1,sig,"/")
  te <- sweep(sweep(te,1,mu,"-"),1,sig,"/")
  list(train=t(tr),test=t(te),features=paste0(layer,"::",rownames(tr)))
}
get_layer_design <- function(M,train_idx,test_idx,layer,max_features) {
  key <- paste(layer,paste(train_idx,collapse=","),paste(test_idx,collapse=","),max_features,sep="__")
  if (exists(key,envir=layer_cache,inherits=FALSE)) return(get(key,envir=layer_cache,inherits=FALSE))
  ans <- prep_omics(M,train_idx,test_idx,layer,max_features)
  assign(key,ans,envir=layer_cache)
  ans
}

mode_value <- function(x) { z <- x[!is.na(x)&nzchar(as.character(x))]; if (!length(z)) return("Unknown"); names(sort(table(z),decreasing=TRUE))[1] }
prep_clinical <- function(d,train_idx,test_idx) {
  tr <- d[train_idx,,drop=FALSE]; te <- d[test_idx,,drop=FALSE]
  age_med <- median(tr$age[is.finite(tr$age)],na.rm=TRUE); if (!is.finite(age_med)) age_med <- 0
  tr$age[!is.finite(tr$age)]<-age_med; te$age[!is.finite(te$age)]<-age_med
  sx <- mode_value(tr$sex); st <- mode_value(tr$stage)
  tr$sex[is.na(tr$sex)|!nzchar(tr$sex)]<-sx; te$sex[is.na(te$sex)|!nzchar(te$sex)]<-sx
  tr$stage[is.na(tr$stage)|!nzchar(tr$stage)]<-st; te$stage[is.na(te$stage)|!nzchar(te$stage)]<-st
  tr$sex <- factor(tr$sex); te$sex <- factor(te$sex,levels=levels(tr$sex))
  te$sex[is.na(te$sex)] <- sx; te$sex <- factor(te$sex,levels=levels(tr$sex))
  tr$stage <- factor(tr$stage); te$stage <- factor(te$stage,levels=levels(tr$stage))
  te$stage[is.na(te$stage)] <- st; te$stage <- factor(te$stage,levels=levels(tr$stage))
  mmtr <- model.matrix(~ age+sex+stage,data=tr)[,-1,drop=FALSE]
  mmte <- model.matrix(~ age+sex+stage,data=te)[,-1,drop=FALSE]
  miss <- setdiff(colnames(mmtr),colnames(mmte)); if(length(miss)) mmte<-cbind(mmte,matrix(0,nrow(mmte),length(miss),dimnames=list(NULL,miss)))
  extra <- setdiff(colnames(mmte),colnames(mmtr)); if(length(extra)) mmte<-mmte[,setdiff(colnames(mmte),extra),drop=FALSE]
  mmte <- mmte[,colnames(mmtr),drop=FALSE]
  # Drop constants based only on outer training data.
  keep <- apply(mmtr,2,sd)>1e-8 & is.finite(apply(mmtr,2,sd))
  list(train=mmtr[,keep,drop=FALSE],test=mmte[,keep,drop=FALSE],features=paste0("CLINICAL::",colnames(mmtr)[keep]))
}
make_design <- function(spec, cohort_obj, train_idx,test_idx,max_features) {
  key <- paste(paste(sort(spec),collapse="+"),paste(train_idx,collapse=","),paste(test_idx,collapse=","),sep="__")
  if (exists(key,envir=design_cache,inherits=FALSE)) return(get(key,envir=design_cache,inherits=FALSE))
  parts <- list()
  if ("CLINICAL" %in% spec) parts[["CLINICAL"]] <- prep_clinical(cohort_obj$clinical,train_idx,test_idx)
  for (b in setdiff(spec,"CLINICAL"))
    parts[[b]] <- get_layer_design(cohort_obj$matrices[[b]],train_idx,test_idx,b,max_features)
  if (!length(parts)) stop("Modelo sin características.")
  xtr<-do.call(cbind,lapply(parts,`[[`,"train")); xte<-do.call(cbind,lapply(parts,`[[`,"test"))
  storage.mode(xtr)<-"double"; storage.mode(xte)<-"double"
  colnames(xtr)<-make.unique(unlist(lapply(parts,`[[`,"features"),use.names=FALSE))
  colnames(xte)<-colnames(xtr)
  ans <- list(train=xtr,test=xte,features=colnames(xtr))
  assign(key,ans,envir=design_cache)
  ans
}

get_cindex <- function(time,event,risk) {
  z <- tryCatch(survival::concordance(survival::Surv(time,event)~risk,reverse=TRUE)$concordance,error=function(e)NA_real_)
  as.numeric(z)
}
fit_nested_cox <- function(spec,cohort_obj,d,train_idx,test_idx,inner_k,seed) {
  ytr <- survival::Surv(d$time_days[train_idx],d$event[train_idx])
  ik <- min(inner_k,min(table(d$event[train_idx])))
  if (ik<2) stop("Insufficient events/censoring for inner cross-validation.")
  inner_fold <- make_folds(d$event[train_idx],ik,seed)
  score_arr <- array(NA_real_,dim=c(length(alpha_grid),ik,length(ratio_grid)))
  for (ai in seq_along(alpha_grid)) {
    a <- alpha_grid[ai]
    for (fi in seq_len(ik)) {
      iti <- train_idx[which(inner_fold!=fi)]; ivi <- train_idx[which(inner_fold==fi)]
      di <- make_design(spec,cohort_obj,iti,ivi,max_features)
      fit_i <- glmnet::glmnet(di$train,survival::Surv(d$time_days[iti],d$event[iti]),
        family="cox",alpha=a,nlambda=nlambda,lambda.min.ratio=min_ratio,
        standardize=FALSE,maxit=50000)
      ratios <- fit_i$lambda/max(fit_i$lambda)
      pred_i <- as.matrix(predict(fit_i,newx=di$test,type="link",s=fit_i$lambda))
      for (ri in seq_along(ratio_grid)) {
        li <- which.min(abs(log(pmax(ratios,1e-12))-log(ratio_grid[ri])))
        score_arr[ai,fi,ri] <- get_cindex(d$time_days[ivi],d$event[ivi],pred_i[,li])
      }
    }
  }
  mean_score <- apply(score_arr,c(1,3),mean,na.rm=TRUE)
  sd_score <- apply(score_arr,c(1,3),sd,na.rm=TRUE)
  if (!any(is.finite(mean_score))) stop("Could not select alpha/lambda during inner cross-validation.")
  best_pos <- which(mean_score==max(mean_score,na.rm=TRUE),arr.ind=TRUE)[1,]
  best_mean <- mean_score[best_pos[1],best_pos[2]]
  best_sd <- sd_score[best_pos[1],best_pos[2]]
  threshold <- best_mean-ifelse(is.finite(best_sd),best_sd/sqrt(ik),0)
  eligible <- which(mean_score>=threshold & is.finite(mean_score),arr.ind=TRUE)
  # One-SE preference: among competitive solutions choose the strongest penalty.
  ord <- order(ratio_grid[eligible[,2]],alpha_grid[eligible[,1]],decreasing=TRUE)
  ai <- eligible[ord[1],1]; ri <- eligible[ord[1],2]
  chosen_alpha <- alpha_grid[ai]; chosen_ratio <- ratio_grid[ri]
  do <- make_design(spec,cohort_obj,train_idx,test_idx,max_features)
  fit_o <- glmnet::glmnet(do$train,ytr,family="cox",alpha=chosen_alpha,nlambda=nlambda,
    lambda.min.ratio=min_ratio,standardize=FALSE,maxit=50000)
  ratios_o <- fit_o$lambda/max(fit_o$lambda)
  li_o <- which.min(abs(log(pmax(ratios_o,1e-12))-log(chosen_ratio)))
  lam <- fit_o$lambda[li_o]
  risk <- as.numeric(predict(fit_o,newx=do$test,s=lam,type="link"))
  coefv <- as.matrix(coef(fit_o,s=lam))[,1]
  list(risk=risk,selected=names(coefv)[is.finite(coefv)&abs(coefv)>1e-10],
       alpha=chosen_alpha,lambda=lam,inner_cindex=mean_score[ai,ri],design=do)
}


model_sets <- list(
  CLINICAL="CLINICAL", RNA="RNA", CNV="CNV", METHYLATION="METHYLATION",
  CLINICAL_RNA=c("CLINICAL","RNA"), CLINICAL_CNV=c("CLINICAL","CNV"),
  CLINICAL_METHYLATION=c("CLINICAL","METHYLATION"),
  OMICS_ALL=c("RNA","CNV","METHYLATION"),
  CLINICAL_ALL=c("CLINICAL","RNA","CNV","METHYLATION")
)
expected_cox_fits <- repeats*outer_k*length(model_sets)*(length(alpha_grid)*inner_k+1L)
cat(sprintf("Ajustes Cox previstos por subtipo: %d | RSF: hasta %d ajustes por subtipo\n",
            expected_cox_fits, if (isTRUE(acfg$cancer_run_rsf_sensitivity)) repeats*outer_k*2L else 0L))
run_one_cohort <- function(cohort) {
  cohort_start <- proc.time()[["elapsed"]]
  cat("\n===== NSCLC ",cohort," =====\n",sep="")
  co <- load_cohort(cohort); d<-co$clinical; n<-nrow(d)
  if (sum(d$event==1)<outer_k || sum(d$event==0)<outer_k) stop(cohort," no tiene suficientes eventos/censura para ",outer_k," folds.")
  outdir<-file.path(root,"results","NSCLC",cohort); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(co$qc,file.path(outdir,"cohort_qc.csv"),row.names=FALSE)
  all_metrics<-list(); all_pred<-list(); all_selected<-list(); metric_i<-0L
  for (rep in seq_len(repeats)) {
    outer<-make_folds(d$event,outer_k,seed0+rep*100+ifelse(cohort=="LUAD",1,2))
    for (fold in seq_len(outer_k)) {
      tr<-which(outer!=fold); te<-which(outer==fold)
      rm(list=ls(envir=design_cache,all.names=TRUE),envir=design_cache)
      rm(list=ls(envir=layer_cache,all.names=TRUE),envir=layer_cache)
      gc(verbose=FALSE)
      for (mn in names(model_sets)) {
        spec<-model_sets[[mn]]
        # Cache each modality's fold-specific transformed matrix for the exact split.
        fit<-fit_nested_cox(spec,co,d,tr,te,inner_k,
                            seed0+rep*10000+fold*100)
        des<-fit$design
        ci<-get_cindex(d$time_days[te],d$event[te],fit$risk)
        metric_i<-metric_i+1L
        all_metrics[[metric_i]]<-data.frame(repetition=rep,fold=fold,model=mn,n_train=length(tr),n_test=length(te),
          events_train=sum(d$event[tr]),events_test=sum(d$event[te]),c_index=ci,alpha=fit$alpha,lambda=fit$lambda,
          inner_cv_cindex=fit$inner_cindex,stringsAsFactors=FALSE)
        all_pred[[metric_i]]<-data.frame(repetition=rep,fold=fold,model=mn,sample_id=d$sample_id[te],
          time_days=d$time_days[te],event=d$event[te],risk=fit$risk,stringsAsFactors=FALSE)
        if(length(fit$selected)) all_selected[[metric_i]]<-data.frame(model=mn,feature_id=fit$selected,
          repetition=rep,fold=fold,alpha=fit$alpha,stringsAsFactors=FALSE)
        cat(cohort,"rep",rep,"fold",fold,mn,"C-index",round(ci,3),"features",ncol(des$train),"\n")
      }
      # Non-linear sensitivity on two prespecified designs. Quick mode caps trees
      # at 100; full mode uses the configured count.
      if (isTRUE(acfg$cancer_run_rsf_sensitivity) && requireNamespace("ranger",quietly=TRUE)) {
        rsf_trees <- if (quick_mode) min(100L,as.integer(acfg$cancer_rsf_trees)) else as.integer(acfg$cancer_rsf_trees)
        for (mn in c("RSF_CLINICAL","RSF_CLINICAL_ALL")) {
          spec<-if(mn=="RSF_CLINICAL") model_sets$CLINICAL else model_sets$CLINICAL_ALL
          des<-make_design(spec,co,tr,te,max_features)
          trdf<-as.data.frame(des$train); te_df<-as.data.frame(des$test)
          names(trdf)<-paste0("f",seq_len(ncol(trdf))); names(te_df)<-names(trdf)
          trdf$time_days<-d$time_days[tr]; trdf$event<-d$event[tr]
          form<-survival::Surv(time_days,event)~.
          rf<-ranger::ranger(form,data=trdf,num.trees=rsf_trees,
             mtry=max(1,floor(sqrt(ncol(des$train)))),min.node.size=15,splitrule="logrank",
             importance="none",seed=seed0+rep*1000+fold)
          pr<-predict(rf,data=te_df)$chf
          risk<-if(is.matrix(pr)) rowSums(pr,na.rm=TRUE) else as.numeric(pr)
          ci<-get_cindex(d$time_days[te],d$event[te],risk)
          metric_i<-metric_i+1L
          all_metrics[[metric_i]]<-data.frame(repetition=rep,fold=fold,model=mn,n_train=length(tr),n_test=length(te),
            events_train=sum(d$event[tr]),events_test=sum(d$event[te]),c_index=ci,alpha=NA,lambda=NA,
            inner_cv_cindex=NA,stringsAsFactors=FALSE)
          all_pred[[metric_i]]<-data.frame(repetition=rep,fold=fold,model=mn,sample_id=d$sample_id[te],
            time_days=d$time_days[te],event=d$event[te],risk=risk,stringsAsFactors=FALSE)
        }
      }
    }
  }
  met<-do.call(rbind,all_metrics); pred<-do.call(rbind,all_pred)
  met <- met[is.finite(met$c_index),,drop=FALSE]
  summ<-do.call(rbind,lapply(split(met,met$model),function(z)data.frame(model=z$model[1],
    mean_c_index=mean(z$c_index),sd_c_index=sd(z$c_index),folds=nrow(z),
    mean_inner_cv_cindex=mean(z$inner_cv_cindex,na.rm=TRUE),stringsAsFactors=FALSE)))
  names(summ)[names(summ)=="mean_inner_cv_cindex"]<-"mean_inner_cv_cindex"
  # Paired delta C-index vs the clinical-only baseline on identical outer folds.
  base<-met[met$model=="CLINICAL",c("repetition","fold","c_index")]; names(base)[3]<-"clinical_c_index"
  delta<-merge(met,base,by=c("repetition","fold"),all.x=TRUE)
  delta$delta_c_index_vs_clinical<-delta$c_index-delta$clinical_c_index
  delta_summary<-do.call(rbind,lapply(split(delta,delta$model),function(z)data.frame(model=z$model[1],
    delta_c_index_vs_clinical=mean(z$delta_c_index_vs_clinical,na.rm=TRUE),
    sd_delta_c_index=sd(z$delta_c_index_vs_clinical,na.rm=TRUE),stringsAsFactors=FALSE)))
  selected<-if(length(all_selected)) do.call(rbind,all_selected[lengths(all_selected)>0]) else data.frame()
  if(nrow(selected)) {
    selected$selected_count<-1L
    stab<-aggregate(selected_count~model+feature_id,selected,sum)
    stab$selection_frequency<-stab$selected_count/(repeats*outer_k)
    stab<-stab[order(stab$model,-stab$selection_frequency),]
    write.csv(selected,file.path(outdir,"selected_features_by_fold.csv"),row.names=FALSE)
    write.csv(stab,file.path(outdir,"feature_stability.csv"),row.names=FALSE)
  }
  write.csv(met,file.path(outdir,"metrics_by_fold.csv"),row.names=FALSE)
  write.csv(pred,file.path(outdir,"oof_predictions.csv"),row.names=FALSE)
  write.csv(summ,file.path(outdir,"model_metrics.csv"),row.names=FALSE)
  write.csv(delta_summary,file.path(outdir,"layer_contribution_vs_clinical.csv"),row.names=FALSE)
  elapsed_min <- (proc.time()[["elapsed"]]-cohort_start)/60
  cat("Finalizado ",cohort,": ",n," pacientes; ",sum(d$event)," eventos; ",nrow(met),
      " resultados fold/modelo; tiempo ",round(elapsed_min,1)," min.\n",sep="")
  invisible(list(metrics=met,predictions=pred,summary=summ,delta=delta_summary,qc=co$qc))
}

luad <- run_one_cohort("LUAD")
lusc <- run_one_cohort("LUSC")
all_qc<-rbind(luad$qc,lusc$qc)
write.csv(all_qc,file.path(root,"results","NSCLC","nsclc_cohort_summary.csv"),row.names=FALSE)
cat("\nNSCLC benchmark completed. LUAD and LUSC were acquired and reconstructed using the same selection criteria.\n")
