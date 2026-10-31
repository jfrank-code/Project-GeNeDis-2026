# Cross-indication gene and GO-term comparison. This is an exploratory transfer
# analysis, not evidence that RA causes cancer or that shared genes have identical effects.
options(stringsAsFactors=FALSE,warn=1)
need<-c("jsonlite","data.table","AnnotationDbi","org.Hs.eg.db","GO.db")
miss<-need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop("Missing packages for cross-indication analysis: ",paste(miss,collapse=", "),". Run src/scripts/R/00_install_dependencies.R")
root<-normalizePath(".",winslash="/",mustWork=TRUE)
out<-file.path(root,"results","CROSS_INDICATION");dir.create(out,recursive=TRUE,showWarnings=FALSE)

map_symbol<-function(ids) {
  x<-sub("\\.[0-9]+$","",as.character(ids)); sym<-x
  ens<-grepl("^ENSG[0-9]+$",x); ent<-grepl("^[0-9]+$",x)
  if(any(ens)) { z<-tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=unique(x[ens]),keytype="ENSEMBL",column="SYMBOL",multiVals="first"),error=function(e)NULL); if(!is.null(z)){ q<-unname(z[match(x,names(z))]);ok<-!is.na(q)&nzchar(q);sym[ok]<-q[ok] } }
  if(any(ent)) { z<-tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=unique(x[ent]),keytype="ENTREZID",column="SYMBOL",multiVals="first"),error=function(e)NULL); if(!is.null(z)){ q<-unname(z[match(x,names(z))]);ok<-!is.na(q)&nzchar(q);sym[ok]<-q[ok] } }
  sym[grepl("^cg[0-9]+$",sym,ignore.case=TRUE)]<-NA_character_
  sym[grepl("^UNMAPPED_",sym)]<-NA_character_
  sym
}
read_stable<-function(path,model_preference=NULL) {
  if(!file.exists(path)) return(data.frame())
  d<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  if(!nrow(d)||!all(c("feature_id","selection_frequency")%in%names(d))) return(data.frame())
  if(!is.null(model_preference)&&"model"%in%names(d)) {
    q<-d[d$model%in%model_preference,,drop=FALSE]
    if(nrow(q)) d<-q
  }
  d$feature_id<-as.character(d$feature_id)
  d$gene_id<-sub("^[^:]+::","",d$feature_id)
  d$gene_symbol<-map_symbol(d$gene_id)
  d<-d[!is.na(d$gene_symbol)&nzchar(d$gene_symbol),,drop=FALSE]
  d
}
# RA prioritizes all-layer selected features, but the output retains their modality.
ra_path<-file.path(root,"results","RA","ra_stable_features.csv")
ra<-read_stable(ra_path,c("ALL","RNA+CHIP","RNA+ATAC"))
if(nrow(ra)) {
  ra<-aggregate(selection_frequency~gene_symbol+modality,ra,function(x)max(x,na.rm=TRUE))
  ra<-ra[order(-ra$selection_frequency,ra$gene_symbol),]
}

nsclc_list<-list()
for(cohort in c("LUAD","LUSC")) {
  p<-file.path(root,"results","NSCLC",cohort,"feature_stability.csv")
  z<-read_stable(p,c("CLINICAL_ALL","OMICS_ALL"))
  if(nrow(z)) {
    z$modality<-sub("::.*$","",z$feature_id)
    z<-z[z$modality%in%c("RNA","CNV"),,drop=FALSE]
    z<-aggregate(selection_frequency~gene_symbol+modality,z,function(x)max(x,na.rm=TRUE))
  }
  nsclc_list[[cohort]]<-z
}

ra_genes<-if(nrow(ra))unique(ra$gene_symbol) else character()
luad_genes<-if(nrow(nsclc_list$LUAD))unique(nsclc_list$LUAD$gene_symbol) else character()
lusc_genes<-if(nrow(nsclc_list$LUSC))unique(nsclc_list$LUSC$gene_symbol) else character()
shared_luad<-intersect(ra_genes,luad_genes); shared_lusc<-intersect(ra_genes,lusc_genes)
all_rows<-list()
add_rows<-function(gset,context,source_df) {
  if(!length(gset)) return(NULL)
  md<-if(nrow(source_df)) aggregate(selection_frequency~gene_symbol,source_df,function(x)max(x,na.rm=TRUE)) else data.frame()
  data.frame(gene_symbol=gset,context=context,
    RA_selection_frequency=if(nrow(ra)) vapply(gset,function(g){v<-ra$selection_frequency[ra$gene_symbol==g];if(length(v))max(v) else NA_real_},numeric(1)) else NA_real_,
    NSCLC_selection_frequency=if(nrow(md)) vapply(gset,function(g){v<-md$selection_frequency[md$gene_symbol==g];if(length(v))max(v) else NA_real_},numeric(1)) else NA_real_,
    stringsAsFactors=FALSE)
}
all_rows[[1]]<-add_rows(shared_luad,"RA_LUAD_shared",nsclc_list$LUAD)
all_rows[[2]]<-add_rows(shared_lusc,"RA_LUSC_shared",nsclc_list$LUSC)
shared<-do.call(rbind,all_rows[!vapply(all_rows,is.null,logical(1))])
if(is.null(shared)||!nrow(shared)) shared<-data.frame(gene_symbol=character(),context=character(),RA_selection_frequency=numeric(),NSCLC_selection_frequency=numeric())
write.csv(shared,file.path(out,"shared_stable_genes.csv"),row.names=FALSE)

# GO Biological Process over-representation. Background is the intersection of measured
# gene symbols between the relevant RA and cancer expression matrices where available.
get_universe<-function(path) {
  if(!file.exists(path)) return(character())
  x<-readRDS(path); if(is.list(x)&&!is.matrix(x)&&!is.data.frame(x)) { q<-x[vapply(x,function(z)is.matrix(z)||is.data.frame(z),logical(1))];if(length(q))x<-q[[1]] }
  map_symbol(rownames(as.matrix(x)))
}
ra_universe<-get_universe(file.path(root,"data","processed","RA","ra_modalities.rds"))
if(!length(ra_universe) && file.exists(file.path(root,"data","processed","RA","ra_modalities.rds"))) {
  o<-readRDS(file.path(root,"data","processed","RA","ra_modalities.rds"));ra_universe<-map_symbol(rownames(o$RNA))
}
get_cancer_universe<-function(cohort) {
  p<-file.path(root,"data","processed",cohort,"RNA.tsv.gz")
  if(!file.exists(p))return(character())
  x<-tryCatch(data.table::fread(p,select=1,data.table=FALSE,showProgress=FALSE)[[1]],error=function(e)character())
  map_symbol(x)
}
# For very large matrices we use a conservative union of measured RNA genes where
# rownames are cheaply available; if not, use all mapped human protein-coding symbols.
run_go<-function(genes,universe,context) {
  genes<-unique(na.omit(genes));universe<-unique(na.omit(universe))
  if(length(genes)<3) return(data.frame())
  if(length(universe)<100) universe<-unique(unlist(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys=AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype="SYMBOL"),keytype="SYMBOL",columns="SYMBOL")$SYMBOL))
  m<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=intersect(universe,genes),keytype="SYMBOL",columns=c("ENTREZID","GO","ONTOLOGY"))
  bg<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=universe,keytype="SYMBOL",columns=c("ENTREZID","GO","ONTOLOGY"))
  bg<-bg[!is.na(bg$GO)&bg$ONTOLOGY=="BP",];m<-m[!is.na(m$GO)&m$ONTOLOGY=="BP",]
  if(!nrow(m)||!nrow(bg)) return(data.frame())
  N<-length(unique(bg$ENTREZID));n<-length(unique(m$ENTREZID));terms<-split(bg$ENTREZID,bg$GO)
  qids<-unique(m$ENTREZID);rows<-lapply(names(terms),function(go){ids<-unique(terms[[go]]);K<-length(ids);k<-length(intersect(qids,ids));if(k<3)return(NULL);data.frame(GO=go,Count=k,UniverseCount=K,pvalue=phyper(k-1,K,N-K,n,lower.tail=FALSE),stringsAsFactors=FALSE)})
  z<-do.call(rbind,rows);if(is.null(z)||!nrow(z))return(data.frame());z$padj<-p.adjust(z$pvalue,"BH");z$Description<-AnnotationDbi::mapIds(GO.db::GO.db,keys=z$GO,keytype="GOID",column="TERM",multiVals="first");z$context<-context;z<-z[order(z$padj,z$pvalue),];head(z,50)
}
# Fallback to a standard human-symbol universe only if exact row names cannot be read.
if(length(ra_universe)<100) ra_universe<-unique(na.omit(ra_genes))
go_rows<-list()
if(length(ra_genes)>=3) go_rows[[length(go_rows)+1]]<-run_go(ra_genes,ra_universe,"RA_stable")
if(length(shared_luad)>=3) go_rows[[length(go_rows)+1]]<-run_go(shared_luad,intersect(ra_universe,get_cancer_universe("LUAD")),"RA_LUAD_shared")
if(length(shared_lusc)>=3) go_rows[[length(go_rows)+1]]<-run_go(shared_lusc,intersect(ra_universe,get_cancer_universe("LUSC")),"RA_LUSC_shared")
go<-if(length(go_rows))do.call(rbind,go_rows) else data.frame()
write.csv(go,file.path(out,"GO_BP_enrichment.csv"),row.names=FALSE)

summary<-data.frame(
  RA_stable_genes=length(ra_genes),LUAD_stable_genes=length(luad_genes),LUSC_stable_genes=length(lusc_genes),
  shared_RA_LUAD=length(shared_luad),shared_RA_LUSC=length(shared_lusc),
  note="Exploratory overlap of stable gene symbols; not causal; genes with unmappable probe IDs excluded.")
write.csv(summary,file.path(out,"cross_indication_summary.csv"),row.names=FALSE)
cat("Cross-indication comparison complete. RA stable:",length(ra_genes),"; shared RA/LUAD:",length(shared_luad),"; shared RA/LUSC:",length(shared_lusc),".\\n")
