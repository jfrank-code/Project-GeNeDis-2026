# MOSAIC-3 | RA paired multi-omics benchmark
# Audited RA multi-omics benchmark implementation
# Task: exploratory RA vs OA classification in donor-matched synovial FLS.
# Models: nested Elastic Net logistic regression + optional linear SVM sensitivity.
# Donor IDs define CV units. All feature filtering, imputation, variance ranking,
# and z-scoring are estimated from each training split only.
# RNA input is the TMM-logCPM matrix produced by 01_prepare_ra_features.R;
# this is label-free normalization performed at preparation time and is disclosed
# in the QC manifest (it is not outcome-guided feature selection).

options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

required_pkgs <- c("glmnet", "jsonlite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

set.seed(20261031)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
infile <- file.path(root, "data", "processed", "RA", "ra_modalities.rds")
manifest_file <- file.path(root, "data", "processed", "RA", "ra_manifest.csv")
out_dir <- file.path(root, "results", "RA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile) || !file.exists(manifest_file)) {
  stop("Missing RA processed inputs. Run src/scripts/R/01_prepare_ra_features.R first.")
}

obj <- readRDS(infile)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("donor", "group") %in% names(manifest))) {
  stop("ra_manifest.csv must contain columns: donor, group.")
}
if (anyDuplicated(manifest$donor)) stop("Duplicate donor IDs in RA manifest.")
manifest$donor <- as.character(manifest$donor)
manifest$group <- factor(toupper(as.character(manifest$group)), levels = c("OA", "RA"))
if (anyNA(manifest$group)) stop("Invalid group labels in RA manifest; expected OA and RA.")

required_blocks <- c("RNA", "ATAC", "H3K27ac", "H3K4me3")
if (!all(required_blocks %in% names(obj))) {
  stop("RDS must contain RNA, ATAC, H3K27ac, and H3K4me3 blocks.")
}

for (b in required_blocks) {
  if (is.null(dim(obj[[b]])) || length(dim(obj[[b]])) != 2L) {
    stop("Omics block is not a 2D matrix: ", b)
  }
  if (is.null(colnames(obj[[b]]))) stop("Missing donor column names in block: ", b)
}

common <- Reduce(intersect, c(list(manifest$donor), lapply(obj[required_blocks], colnames)))
if (length(common) < 12L) stop("Fewer than 12 donors have all four paired modalities.")
manifest <- manifest[match(common, manifest$donor), , drop = FALSE]
if (anyNA(manifest$donor) || anyDuplicated(manifest$donor)) stop("Donor alignment failed.")
if (length(unique(manifest$group)) != 2L) stop("Both RA and OA groups are required.")
obj <- lapply(obj[required_blocks], function(m) m[, common, drop = FALSE])
y <- as.integer(manifest$group == "RA")

cfg_path <- file.path(root, "config", "mosaic3_config.json")
if (!file.exists(cfg_path)) stop("Missing config/mosaic3_config.json")
cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
a <- cfg$analysis
repeats <- as.integer(a$ra_outer_repeats %||% 2L)
outer_k <- as.integer(a$ra_outer_folds %||% 5L)
inner_k <- as.integer(a$ra_inner_folds %||% 3L)
max_features <- as.integer(a$ra_max_features_per_layer %||% 120L)
nlambda <- as.integer(a$ra_nlambda %||% 10L)
min_ratio <- as.numeric(a$ra_lambda_min_ratio %||% 0.05)
run_svm <- isTRUE(a$ra_run_svm_sensitivity) && requireNamespace("e1071", quietly = TRUE)
alpha <- 0.5

if (any(!is.finite(c(repeats, outer_k, inner_k, max_features, nlambda, min_ratio))) ||
    repeats < 1L || outer_k < 2L || inner_k < 2L || max_features < 1L || nlambda < 3L ||
    min_ratio <= 0 || min_ratio >= 1) {
  stop("Invalid RA analysis settings in mosaic3_config.json.")
}
if (min(table(y)) < outer_k) stop("Not enough donors per class for requested outer folds.")

cat("RA paired multi-omics benchmark\n")
cat("Donors:", length(common), "| RA:", sum(y == 1L), "| OA:", sum(y == 0L), "\n")
cat("Repeats:", repeats, "| outer folds:", outer_k, "| inner folds:", inner_k,
    "| max features/layer:", max_features, "| lambdas:", nlambda,
    "| SVM sensitivity:", run_svm, "\n")


make_folds <- function(y_vec, k, seed) {
  set.seed(seed)
  f <- integer(length(y_vec))
  for (cl in sort(unique(y_vec))) {
    ii <- which(y_vec == cl)
    ii <- sample(ii)
    f[ii] <- rep(seq_len(k), length.out = length(ii))
  }
  if (any(f == 0L)) stop("Could not construct stratified folds.")
  f
}

auc_rank <- function(y_true, score) {
  ok <- is.finite(score) & !is.na(y_true)
  y_true <- as.integer(y_true[ok])
  score <- as.numeric(score[ok])
  n1 <- sum(y_true == 1L)
  n0 <- sum(y_true == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y_true == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_metrics <- function(y_true, score, threshold = 0.5, is_probability = TRUE) {
  ok <- is.finite(score) & !is.na(y_true)
  y_true <- as.integer(y_true[ok])
  score <- as.numeric(score[ok])
  if (!length(y_true)) {
    return(c(AUC = NA_real_, BalancedAccuracy = NA_real_, Sensitivity = NA_real_,
             Specificity = NA_real_, MCC = NA_real_, Brier = NA_real_))
  }
  pred <- as.integer(score >= threshold)
  tp <- sum(pred == 1L & y_true == 1L)
  tn <- sum(pred == 0L & y_true == 0L)
  fp <- sum(pred == 1L & y_true == 0L)
  fn <- sum(pred == 0L & y_true == 1L)
  sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (is.finite(mcc_den) && mcc_den > 0) (tp * tn - fp * fn) / mcc_den else NA_real_
  ba <- if (is.finite(sens) && is.finite(spec)) mean(c(sens, spec)) else NA_real_
  brier <- if (is_probability) mean((score - y_true)^2) else NA_real_
  c(AUC = auc_rank(y_true, score), BalancedAccuracy = ba,
    Sensitivity = sens, Specificity = spec, MCC = mcc, Brier = brier)
}

# Track model warnings rather than flooding the console or silently discarding them.
warning_log <- data.frame(context = character(), warning = character(), stringsAsFactors = FALSE)
record_warning <- function(context, message) {
  warning_log <<- rbind(warning_log, data.frame(context = as.character(context),
                                               warning = as.character(message),
                                               stringsAsFactors = FALSE))
}

# Fold-specific preprocessing: all medians, variance ranking and scaling are
# estimated from training samples only. Input matrices are features x donors.
prep_block <- function(M, train_idx, test_idx, layer, max_features) {
  M <- as.matrix(M)
  suppressWarnings(storage.mode(M) <- "numeric")
  if (is.null(rownames(M))) rownames(M) <- paste0(layer, "_feature_", seq_len(nrow(M)))
  tr <- M[, train_idx, drop = FALSE]
  te <- M[, test_idx, drop = FALSE]

  if (layer != "RNA") {
    tr <- log1p(pmax(tr, 0))
    te <- log1p(pmax(te, 0))
  }

  if (nrow(tr) == 0L) stop("No features available for layer: ", layer)
  med <- apply(tr, 1L, median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(nrow(tr))) {
    tr[j, !is.finite(tr[j, ])] <- med[j]
    te[j, !is.finite(te[j, ])] <- med[j]
  }

  vars <- apply(tr, 1L, var)
  vars[!is.finite(vars)] <- 0
  ord <- order(vars, decreasing = TRUE)
  ord <- ord[vars[ord] > 1e-10]
  if (!length(ord)) stop("No variable features in layer: ", layer)
  ord <- head(ord, max_features)
  tr <- tr[ord, , drop = FALSE]
  te <- te[ord, , drop = FALSE]

  mu <- rowMeans(tr)
  sig <- apply(tr, 1L, sd)
  keep <- is.finite(sig) & sig > 1e-8
  tr <- tr[keep, , drop = FALSE]
  te <- te[keep, , drop = FALSE]
  mu <- mu[keep]
  sig <- sig[keep]
  if (!nrow(tr)) stop("All selected features had zero variance in layer: ", layer)

  tr <- sweep(sweep(tr, 1L, mu, "-"), 1L, sig, "/")
  te <- sweep(sweep(te, 1L, mu, "-"), 1L, sig, "/")
  colnames(tr) <- colnames(M)[train_idx]
  colnames(te) <- colnames(M)[test_idx]
  list(x_train = t(tr), x_test = t(te), features = paste0(layer, "::", rownames(tr)))
}

layer_cache <- new.env(parent = emptyenv())
get_prepped_block <- function(block, train_idx, test_idx) {
  key <- paste(block, paste(train_idx, collapse = ","), paste(test_idx, collapse = ","), sep = "__")
  if (exists(key, envir = layer_cache, inherits = FALSE)) {
    return(get(key, envir = layer_cache, inherits = FALSE))
  }
  value <- prep_block(obj[[block]], train_idx, test_idx, block, max_features)
  assign(key, value, envir = layer_cache)
  value
}

make_design <- function(block_names, train_idx, test_idx) {
  parts <- lapply(block_names, function(b) get_prepped_block(b, train_idx, test_idx))
  x_train <- do.call(cbind, lapply(parts, `[[`, "x_train"))
  x_test <- do.call(cbind, lapply(parts, `[[`, "x_test"))
  if (is.null(dim(x_train)) || ncol(x_train) == 0L) stop("Empty training design matrix.")
  # Preserve the modality namespace in the actual glmnet column names.
  # This prevents same-symbol features from different omics layers being conflated.
  feature_names <- make.unique(unlist(lapply(parts, `[[`, "features"), use.names = FALSE))
  if (length(feature_names) != ncol(x_train) || ncol(x_train) != ncol(x_test)) {
    stop("Feature-name/design-matrix dimension mismatch.")
  }
  colnames(x_train) <- feature_names
  colnames(x_test) <- feature_names
  list(x_train = x_train, x_test = x_test, features = feature_names)
}

ratio_grid <- exp(seq(log(1), log(min_ratio), length.out = nlambda))
fit_en <- function(x_train, y_train, x_test, context) {
  captured <- character()
  fit <- withCallingHandlers(
    glmnet::glmnet(x = x_train, y = y_train, family = "binomial", alpha = alpha,
                   nlambda = nlambda, lambda.min.ratio = min_ratio,
                   standardize = FALSE, maxit = 50000),
    warning = function(w) {
      captured <<- c(captured, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(captured)) {
    for (w in unique(captured)) record_warning(context, w)
  }
  pred <- as.matrix(predict(fit, newx = x_test, type = "response"))
  actual_ratio <- fit$lambda / max(fit$lambda)
  list(fit = fit, pred = pred, ratios = actual_ratio)
}

sets <- list(
  RNA = "RNA",
  ATAC = "ATAC",
  H3K27ac = "H3K27ac",
  H3K4me3 = "H3K4me3",
  `RNA+ATAC` = c("RNA", "ATAC"),
  `RNA+CHIP` = c("RNA", "H3K27ac", "H3K4me3"),
  `ATAC+CHIP` = c("ATAC", "H3K27ac", "H3K4me3"),
  ALL = required_blocks
)

all_pred <- list()
all_metrics <- list()
all_selected <- list()
feature_sets <- list()
completed_model_folds <- list()
completed_fold_counter <- 0L
pred_counter <- 0L
metric_counter <- 0L
selected_counter <- 0L
feature_set_counter <- 0L

for (rep_id in seq_len(repeats)) {
  outer_fold <- make_folds(y, outer_k, 20261031 + rep_id)
  rep_pred <- list()

  for (fold_id in seq_len(outer_k)) {
    rm(list = ls(envir = layer_cache, all.names = TRUE), envir = layer_cache)
    train_idx <- which(outer_fold != fold_id)
    test_idx <- which(outer_fold == fold_id)
    y_train <- y[train_idx]
    inner_fold <- make_folds(y_train, inner_k, 20261031 + rep_id * 100L + fold_id)

    for (model_name in names(sets)) {
      blocks <- sets[[model_name]]
      inner_scores <- matrix(NA_real_, nrow = inner_k, ncol = nlambda)

      for (inner_id in seq_len(inner_k)) {
        inner_train_local <- which(inner_fold != inner_id)
        inner_valid_local <- which(inner_fold == inner_id)
        inner_train_idx <- train_idx[inner_train_local]
        inner_valid_idx <- train_idx[inner_valid_local]

        # A validation fold with one class cannot estimate AUC; skip that fold.
        if (length(unique(y[inner_valid_idx])) < 2L || length(unique(y[inner_train_idx])) < 2L) {
          record_warning(paste("inner", rep_id, fold_id, model_name),
                         "Inner split had only one class; skipped for lambda selection.")
          next
        }

        design_inner <- tryCatch(
          make_design(blocks, inner_train_idx, inner_valid_idx),
          error = function(e) {
            record_warning(paste("inner", rep_id, fold_id, model_name), conditionMessage(e))
            NULL
          }
        )
        if (is.null(design_inner)) next

        fit_inner <- tryCatch(
          fit_en(design_inner$x_train, y[inner_train_idx], design_inner$x_test,
                 paste("inner", rep_id, fold_id, model_name, inner_id)),
          error = function(e) {
            record_warning(paste("inner", rep_id, fold_id, model_name), conditionMessage(e))
            NULL
          }
        )
        if (is.null(fit_inner)) next

        vals <- vapply(seq_len(ncol(fit_inner$pred)), function(j) {
          auc_rank(y[inner_valid_idx], fit_inner$pred[, j])
        }, numeric(1))
        ratios <- fit_inner$ratios
        valid <- is.finite(vals) & is.finite(ratios)
        if (sum(valid) >= 2L) {
          ord <- order(log(ratios[valid]))
          inner_scores[inner_id, ] <- approx(
            x = log(ratios[valid])[ord], y = vals[valid][ord],
            xout = log(ratio_grid), rule = 2, ties = mean
          )$y
        }
      }

      n_valid_inner <- colSums(is.finite(inner_scores))
      mean_auc <- apply(inner_scores, 2L, function(z) {
        if (all(!is.finite(z))) NA_real_ else mean(z, na.rm = TRUE)
      })
      sd_auc <- apply(inner_scores, 2L, function(z) {
        if (sum(is.finite(z)) < 2L) NA_real_ else sd(z, na.rm = TRUE)
      })
      if (!any(is.finite(mean_auc))) {
        record_warning(paste("outer", rep_id, fold_id, model_name),
                       "No valid inner AUC values; model/fold skipped.")
        next
      }

      best_idx <- which.max(mean_auc)
      se <- if (is.finite(sd_auc[best_idx])) {
        sd_auc[best_idx] / sqrt(max(1L, n_valid_inner[best_idx]))
      } else 0
      threshold <- mean_auc[best_idx] - se
      eligible <- which(is.finite(mean_auc) & mean_auc >= threshold)
      selected_ratio <- max(ratio_grid[eligible]) # stronger penalty within one-SE set

      design_outer <- tryCatch(
        make_design(blocks, train_idx, test_idx),
        error = function(e) {
          record_warning(paste("outer", rep_id, fold_id, model_name), conditionMessage(e))
          NULL
        }
      )
      if (is.null(design_outer)) next

      fit_outer <- tryCatch(
        fit_en(design_outer$x_train, y_train, design_outer$x_test,
               paste("outer", rep_id, fold_id, model_name)),
        error = function(e) {
          record_warning(paste("outer", rep_id, fold_id, model_name), conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit_outer)) next

      lambda_idx <- which.min(abs(log(pmax(fit_outer$ratios, 1e-12)) - log(selected_ratio)))
      prob <- as.numeric(fit_outer$pred[, lambda_idx])
      if (length(prob) != length(test_idx)) {
        record_warning(paste("outer", rep_id, fold_id, model_name),
                       "Prediction length mismatch; model/fold skipped.")
        next
      }

      rep_pred[[length(rep_pred) + 1L]] <- data.frame(
        repetition = rep_id, fold = fold_id, model = model_name,
        donor = manifest$donor[test_idx], group = as.character(manifest$group[test_idx]),
        truth = y[test_idx], probability = prob, stringsAsFactors = FALSE
      )
      completed_fold_counter <- completed_fold_counter + 1L
      completed_model_folds[[completed_fold_counter]] <- data.frame(
        model = model_name, repetition = rep_id, fold = fold_id,
        stringsAsFactors = FALSE
      )

      co <- as.matrix(stats::coef(fit_outer$fit, s = fit_outer$fit$lambda[lambda_idx]))[, 1L]
      nz <- which(is.finite(co) & abs(co) > 1e-10 & names(co) != "(Intercept)")
      selected_names <- if (length(nz)) names(co)[nz] else character()
      feature_set_counter <- feature_set_counter + 1L
      feature_sets[[feature_set_counter]] <- list(
        model = model_name, repetition = rep_id, fold = fold_id,
        feature_id = selected_names
      )
      if (length(nz)) {
        selected_counter <- selected_counter + 1L
        all_selected[[selected_counter]] <- data.frame(
          model = model_name, feature_id = names(co)[nz], coefficient = as.numeric(co[nz]),
          repetition = rep_id, fold = fold_id, stringsAsFactors = FALSE
        )
      }
      cat("RA rep", rep_id, "fold", fold_id, "model", model_name,
          "features", ncol(design_outer$x_train), "\n")
    }
  }

  if (!length(rep_pred)) stop("No RA models produced predictions in repetition ", rep_id)
  rp <- do.call(rbind, rep_pred)
  pred_counter <- pred_counter + 1L
  all_pred[[pred_counter]] <- rp

  # Pooled OOF metrics: each donor is tested once per model in each repetition.
  for (model_name in unique(rp$model)) {
    z <- rp[rp$model == model_name, , drop = FALSE]
    mm <- calc_metrics(z$truth, z$probability)
    metric_counter <- metric_counter + 1L
    all_metrics[[metric_counter]] <- data.frame(
      repetition = rep_id, fold = "pooled_OOF", model = model_name,
      t(mm), n_test = nrow(z), check.names = FALSE
    )
  }
  # Checkpoint after each completed outer repetition to protect work if interrupted.
  utils::write.csv(do.call(rbind, all_pred), file.path(out_dir, "ra_oof_predictions_partial.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, all_metrics), file.path(out_dir, "ra_metrics_partial.csv"), row.names = FALSE)
  cat("Checkpoint saved after RA repetition", rep_id, "of", repeats, "\n")
}

# Fixed linear SVM sensitivity analysis. The orientation is learned ONLY from
# training-set decision values and training labels, never from held-out labels.
if (run_svm) {
  for (rep_id in seq_len(repeats)) {
    outer_fold <- make_folds(y, outer_k, 20261031 + rep_id)
    svm_rep <- list()

    svm_specs <- list(RNA_SVM = "RNA", ALL_SVM = required_blocks)
    for (fold_id in seq_len(outer_k)) {
      rm(list = ls(envir = layer_cache, all.names = TRUE), envir = layer_cache)
      train_idx <- which(outer_fold != fold_id)
      test_idx <- which(outer_fold == fold_id)

      for (nm in names(svm_specs)) {
        blocks <- svm_specs[[nm]]
        des <- tryCatch(make_design(blocks, train_idx, test_idx), error = function(e) {
          record_warning(paste("SVM", rep_id, fold_id, nm), conditionMessage(e))
          NULL
        })
        if (is.null(des)) next

        fit <- tryCatch(
          e1071::svm(
            x = des$x_train,
            y = factor(y[train_idx], levels = c(0, 1)),
            kernel = "linear", cost = 0.5, scale = FALSE
          ),
          error = function(e) {
            record_warning(paste("SVM", rep_id, fold_id, nm), conditionMessage(e))
            NULL
          }
        )
        if (is.null(fit)) next

        get_decision <- function(model, newdata, expected_n) {
          pred_obj <- predict(model, newdata, decision.values = TRUE)
          dv <- attr(pred_obj, "decision.values")
          if (is.null(dv)) return(NULL)
          dv <- as.matrix(dv)
          if (nrow(dv) != expected_n || ncol(dv) != 1L) return(NULL)
          score <- as.numeric(dv[, 1L])
          if (length(score) != expected_n || any(!is.finite(score))) return(NULL)
          score
        }

        # Correctly obtain a training score vector of length n_train.
        train_score <- tryCatch(get_decision(fit, des$x_train, length(train_idx)), error = function(e) NULL)
        test_score <- tryCatch(get_decision(fit, des$x_test, length(test_idx)), error = function(e) NULL)
        if (is.null(train_score) || is.null(test_score)) {
          record_warning(paste("SVM", rep_id, fold_id, nm),
                         "Could not obtain finite binary decision values; SVM fold skipped.")
          next
        }

        mean_ra <- mean(train_score[y[train_idx] == 1L])
        mean_oa <- mean(train_score[y[train_idx] == 0L])
        if (!is.finite(mean_ra) || !is.finite(mean_oa) || abs(mean_ra - mean_oa) < 1e-12) {
          record_warning(paste("SVM", rep_id, fold_id, nm),
                         "Decision-score orientation was ambiguous; SVM fold skipped.")
          next
        }
        if (mean_ra < mean_oa) test_score <- -test_score

        svm_rep[[length(svm_rep) + 1L]] <- data.frame(
          repetition = rep_id, fold = fold_id, model = nm,
          donor = manifest$donor[test_idx], group = as.character(manifest$group[test_idx]),
          truth = y[test_idx], probability = test_score, stringsAsFactors = FALSE
        )
        cat("RA SVM rep", rep_id, "fold", fold_id, nm,
            "features", ncol(des$x_train), "\n")
      }
    }

    if (length(svm_rep)) {
      zrep <- do.call(rbind, svm_rep)
      for (nm in unique(zrep$model)) {
        z <- zrep[zrep$model == nm, , drop = FALSE]
        all_pred[[length(all_pred) + 1L]] <- z
        mm <- calc_metrics(z$truth, z$probability, threshold = 0, is_probability = FALSE)
        all_metrics[[length(all_metrics) + 1L]] <- data.frame(
          repetition = rep_id, fold = "pooled_OOF", model = nm,
          t(mm), n_test = nrow(z), check.names = FALSE
        )
      }
    }
  }
} else if (isTRUE(a$ra_run_svm_sensitivity)) {
  record_warning("SVM", "e1071 is not installed; SVM sensitivity analysis was skipped.")
}

if (!length(all_pred) || !length(all_metrics)) stop("No RA predictions or metrics were generated.")
pred <- do.call(rbind, all_pred)
met <- do.call(rbind, all_metrics)

metric_names <- c("AUC", "BalancedAccuracy", "Sensitivity", "Specificity", "MCC", "Brier")
metric_mean <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
metric_sd <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)

summary <- do.call(rbind, lapply(split(met, met$model), function(z) {
  row <- data.frame(model = z$model[1], n_repetitions = nrow(z), stringsAsFactors = FALSE)
  for (nm in metric_names) {
    safe_nm <- tolower(gsub("([a-z])([A-Z])", "\\1_\\2", nm))
    row[[paste0("mean_", safe_nm)]] <- metric_mean(z[[nm]])
    row[[paste0("sd_", safe_nm)]] <- metric_sd(z[[nm]])
  }
  row
}))

repeat_summary <- do.call(rbind, lapply(split(met, paste(met$repetition, met$model, sep = "__")), function(z) {
  row <- data.frame(repetition = z$repetition[1], model = z$model[1], stringsAsFactors = FALSE)
  for (nm in metric_names) {
    safe_nm <- tolower(gsub("([a-z])([A-Z])", "\\1_\\2", nm))
    row[[paste0("mean_", safe_nm)]] <- metric_mean(z[[nm]])
  }
  row
}))

# Incremental contribution relative to RNA-only, paired within the same repeat.
base_metrics <- met[met$model == "RNA", c("repetition", "AUC", "BalancedAccuracy"), drop = FALSE]
names(base_metrics) <- c("repetition", "rna_auc", "rna_balanced_accuracy")
contrib <- merge(met, base_metrics, by = "repetition", all.x = TRUE)
contrib$delta_auc_vs_RNA <- contrib$AUC - contrib$rna_auc
contrib$delta_balanced_accuracy_vs_RNA <- contrib$BalancedAccuracy - contrib$rna_balanced_accuracy
contrib_out <- do.call(rbind, lapply(split(contrib, contrib$model), function(z) data.frame(
  model = z$model[1],
  delta_auc_vs_RNA = metric_mean(z$delta_auc_vs_RNA),
  sd_delta_auc_vs_RNA = metric_sd(z$delta_auc_vs_RNA),
  delta_balanced_accuracy_vs_RNA = metric_mean(z$delta_balanced_accuracy_vs_RNA),
  sd_delta_balanced_accuracy_vs_RNA = metric_sd(z$delta_balanced_accuracy_vs_RNA),
  stringsAsFactors = FALSE
)))

# Biomarker stability: selection frequency and pairwise Jaccard across outer folds.
selected_out <- data.frame()
stability_out <- data.frame()
jaccard_out <- data.frame()
if (length(all_selected)) {
  selected_out <- do.call(rbind, all_selected)
  # One feature can only contribute once per model/repetition/fold. Keep distinct
  # modality-qualified IDs separate, but remove duplicate rows from the same layer.
  selected_out <- selected_out[!duplicated(selected_out[c("model", "feature_id", "repetition", "fold")]), , drop = FALSE]
  selected_out$modality <- sub("::.*$", "", selected_out$feature_id)
  selected_out$gene_or_feature <- sub("^[^:]+::", "", selected_out$feature_id)
  selected_out$selected_count <- 1L
  stability_out <- aggregate(selected_count ~ model + feature_id + modality + gene_or_feature,
                             selected_out, sum)
  completed_df <- unique(do.call(rbind, completed_model_folds))
  denom <- aggregate(fold ~ model, completed_df, length)
  names(denom)[2] <- "completed_outer_folds"
  stability_out <- merge(stability_out, denom, by = "model", all.x = TRUE)
  stability_out$selection_frequency <- stability_out$selected_count / stability_out$completed_outer_folds
  stability_out <- stability_out[order(stability_out$model, -stability_out$selection_frequency), ]
}

if (length(feature_sets)) {
  model_names <- unique(vapply(feature_sets, function(x) x$model, character(1)))
  for (nm in model_names) {
    items <- Filter(function(x) identical(x$model, nm), feature_sets)
    feature_lists <- lapply(items, function(x) unique(x$feature_id[!is.na(x$feature_id) & nzchar(x$feature_id)]))
    if (length(feature_lists) >= 2L) {
      pairs <- utils::combn(seq_along(feature_lists), 2L)
      vals <- apply(pairs, 2L, function(ix) {
        aa <- feature_lists[[ix[1L]]]
        bb <- feature_lists[[ix[2L]]]
        uu <- union(aa, bb)
        if (!length(uu)) NA_real_ else length(intersect(aa, bb)) / length(uu)
      })
      jaccard_out <- rbind(jaccard_out, data.frame(
        model = nm, comparisons = sum(is.finite(vals)),
        mean_jaccard = metric_mean(vals), sd_jaccard = metric_sd(vals),
        stringsAsFactors = FALSE
      ))
    }
  }
}

utils::write.csv(contrib_out, file.path(out_dir, "ra_layer_contribution.csv"), row.names = FALSE)
utils::write.csv(pred, file.path(out_dir, "ra_oof_predictions.csv"), row.names = FALSE)
utils::write.csv(met, file.path(out_dir, "ra_metrics_by_oof_repeat.csv"), row.names = FALSE)
utils::write.csv(repeat_summary, file.path(out_dir, "ra_metrics_by_repeat.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "ra_model_metrics.csv"), row.names = FALSE)
if (nrow(selected_out)) utils::write.csv(selected_out, file.path(out_dir, "ra_selected_features_by_fold.csv"), row.names = FALSE)
if (nrow(stability_out)) utils::write.csv(stability_out, file.path(out_dir, "ra_stable_features.csv"), row.names = FALSE)
if (nrow(jaccard_out)) utils::write.csv(jaccard_out, file.path(out_dir, "ra_stability_jaccard.csv"), row.names = FALSE)
if (nrow(warning_log)) utils::write.csv(warning_log, file.path(out_dir, "ra_model_warnings.csv"), row.names = FALSE)

qc_manifest <- data.frame(
  donors = length(common), RA = sum(y == 1L), OA = sum(y == 0L),
  outer_repeats = repeats, outer_folds = outer_k, inner_folds = inner_k,
  max_features_per_layer = max_features, nlambda = nlambda, alpha = alpha,
  svm_sensitivity_run = run_svm,
  preprocessing = paste(
    "RNA TMM-logCPM prepared once without labels; fold-specific median imputation,",
    "variance ranking, feature filtering, and z-scoring; ATAC/ChIP log1p transform;",
    "all supervised feature selection and lambda tuning restricted to training folds."
  ),
  caveat = paste(
    "Small donor cohort: inner-training class counts can be below glmnet's recommended",
    "threshold of 8; warnings are recorded. Treat all RA metrics as exploratory."
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(qc_manifest, file.path(out_dir, "ra_qc_model_manifest.csv"), row.names = FALSE)

# Remove transient checkpoints only after final outputs were written successfully.
unlink(file.path(out_dir, c("ra_oof_predictions_partial.csv", "ra_metrics_partial.csv")))

cat("\nRA benchmark completed. Donors:", length(common),
    "| predictions:", nrow(pred), "| metric rows:", nrow(met), "\n")
cat("Elastic Net modality configurations:", length(sets),
    "| SVM sensitivity:", run_svm, "\n")
cat("Warnings recorded:", nrow(warning_log), "\n")
cat("Outputs written to:", out_dir, "\n")
