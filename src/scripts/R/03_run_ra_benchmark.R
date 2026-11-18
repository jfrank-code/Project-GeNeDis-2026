# MOSAIC-3 | RA/OA donor-level nested multi-omics classification benchmark
# Leakage-conscious revision for the small paired synovial cohort.
#
# Key safeguards and changes:
# - Donor is the CV unit; no donor is split across train/test within a fold.
# - RNA is stored as raw counts. TMM reference factors and logCPM are estimated
#   within each training split; each held-out donor is normalized relative to a
#   fixed training reference donor (no cohort-wide TMM before CV).
# - ATAC/H3K27ac/H3K4me3 features use fixed promoter universe and donor-local
#   signal normalization performed in 01_prepare_ra_features.R.
# - Expression/prevalence filters, imputation, feature ranking, and scaling are
#   estimated using training donors only (including inner CV training partitions).
# - Feature count is capped by model size: at most 30 total features per model,
#   and at most 20 for a single-layer model.
# - Lambda selection uses pooled inner out-of-fold AUC across inner validation
#   donors, then chooses the strongest penalty within 0.02 AUC of the best score.
# - Outer-test labels are used only once for final evaluation, never for tuning,
#   preprocessing, orientation, or threshold selection.
#
# All metrics remain exploratory because the cohort has only 20 donors.

options(stringsAsFactors = FALSE, warn = 1, timeout = 3600)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

required_pkgs <- c("glmnet", "jsonlite", "edgeR", "matrixStats")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

set.seed(20261117)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
infile <- file.path(root, "data", "processed", "RA", "ra_modalities.rds")
manifest_file <- file.path(root, "data", "processed", "RA", "ra_manifest.csv")
out_dir <- file.path(root, "results", "RA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile) || !file.exists(manifest_file)) {
  stop("Missing RA processed inputs. Run 01_prepare_ra_features.R first.")
}
obj <- readRDS(infile)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("donor", "group") %in% names(manifest))) {
  stop("ra_manifest.csv must contain donor and group columns.")
}
if (anyDuplicated(manifest$donor)) stop("Duplicate donor IDs in RA manifest.")
manifest$donor <- as.character(manifest$donor)
manifest$group <- factor(toupper(as.character(manifest$group)), levels = c("OA", "RA"))
if (anyNA(manifest$group)) stop("Invalid group labels; expected RA and OA.")

required_blocks <- c("RNA", "ATAC", "H3K27ac", "H3K4me3")
if (!all(required_blocks %in% names(obj))) stop("RDS must contain RNA, ATAC, H3K27ac and H3K4me3.")
for (b in required_blocks) {
  if (is.null(dim(obj[[b]])) || length(dim(obj[[b]])) != 2L) stop("Invalid matrix block: ", b)
  if (is.null(colnames(obj[[b]]))) stop("Missing donor column names in block: ", b)
  if (anyDuplicated(rownames(obj[[b]]))) stop("Duplicated feature IDs in block: ", b)
}

common <- Reduce(intersect, c(list(manifest$donor), lapply(obj[required_blocks], colnames)))
common <- sort(common)
if (length(common) < 12L) stop("Fewer than 12 donors have all four paired modalities.")
manifest <- manifest[match(common, manifest$donor), , drop = FALSE]
if (anyNA(manifest$donor) || anyDuplicated(manifest$donor)) stop("Donor alignment failed.")
if (length(unique(manifest$group)) != 2L) stop("Both RA and OA groups are required.")
obj <- lapply(obj[required_blocks], function(m) m[, common, drop = FALSE])
# Guard against accidentally reusing the old processed RDS, which stored
# cohort-wide TMM-logCPM values. The revised benchmark requires raw integer RNA counts.
rna_values <- as.matrix(obj$RNA)
if (any(!is.finite(rna_values)) || any(rna_values < 0) ||
    any(abs(rna_values - round(rna_values)) > 1e-6)) {
  stop("The RA RDS does not contain raw integer RNA counts. Re-run the revised ",
       "01_prepare_ra_features.R before running this benchmark.")
}
y <- as.integer(manifest$group == "RA")
if (min(table(y)) < 2L) stop("At least two donors per class are required.")

cfg_path <- file.path(root, "config", "mosaic3_config.json")
if (!file.exists(cfg_path)) stop("Missing config/mosaic3_config.json")
cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
a <- cfg$analysis
configured_repeats <- as.integer(a$ra_outer_repeats %||% 10L)
# Allow up to ten outer repetitions to examine sensitivity to donor partitions.
# Repeated CV reuses the same 20 donors; it does NOT add independent samples
# and must NOT be interpreted as external validation.
if (is.na(configured_repeats) || configured_repeats < 1L) {
  stop("analysis.ra_outer_repeats must be a positive integer in config/mosaic3_config.json.")
}
repeats <- min(max(3L, configured_repeats), 20L)
outer_k <- min(as.integer(a$ra_outer_folds %||% 5L), 5L)
inner_k_cfg <- as.integer(a$ra_inner_folds %||% 3L)
# Keep the inner split count feasible for tiny training folds.
inner_k <- min(inner_k_cfg, 3L)
configured_max <- as.integer(a$ra_max_features_per_layer %||% 20L)
max_features_per_layer <- min(configured_max, 20L)
nlambda <- min(max(as.integer(a$ra_nlambda %||% 15L), 15L), 20L)
min_ratio <- as.numeric(a$ra_lambda_min_ratio %||% 0.03)
run_svm <- isTRUE(a$ra_run_svm_sensitivity) && requireNamespace("e1071", quietly = TRUE)
alpha <- 0.5
max_total_features <- 30L
inner_auc_tolerance <- 0.02

if (any(!is.finite(c(repeats, outer_k, inner_k, max_features_per_layer, nlambda, min_ratio))) ||
    repeats < 1L || outer_k < 2L || inner_k < 2L || max_features_per_layer < 5L ||
    nlambda < 5L || min_ratio <= 0 || min_ratio >= 1) {
  stop("Invalid RA settings in config/mosaic3_config.json.")
}
if (min(table(y)) < outer_k) stop("Not enough donors per class for requested outer folds.")

cat("RA paired multi-omics benchmark | leakage-conscious revision\n")
cat("Donors:", length(common), "| RA:", sum(y == 1L), "| OA:", sum(y == 0L), "\n")
cat("Outer repeats:", repeats, "| outer folds:", outer_k, "| inner folds:", inner_k,
    "| max total features/model:", max_total_features, "| max/layer:", max_features_per_layer,
    "| lambda grid:", nlambda, "| SVM sensitivity:", run_svm, "\n")

make_folds <- function(y_vec, k, seed) {
  if (length(y_vec) < k || min(table(y_vec)) < k) stop("Not enough observations per class for ", k, " stratified folds.")
  set.seed(seed)
  fold <- integer(length(y_vec))
  for (cl in sort(unique(y_vec))) {
    ix <- sample(which(y_vec == cl))
    fold[ix] <- rep(seq_len(k), length.out = length(ix))
  }
  if (any(fold == 0L)) stop("Could not construct stratified folds.")
  fold
}

auc_rank <- function(y_true, score) {
  ok <- is.finite(score) & !is.na(y_true)
  y_true <- as.integer(y_true[ok]); score <- as.numeric(score[ok])
  n1 <- sum(y_true == 1L); n0 <- sum(y_true == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y_true == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_metrics <- function(y_true, score, threshold = 0.5, is_probability = TRUE) {
  ok <- is.finite(score) & !is.na(y_true)
  y_true <- as.integer(y_true[ok]); score <- as.numeric(score[ok])
  if (!length(y_true)) return(c(AUC = NA_real_, BalancedAccuracy = NA_real_,
    Sensitivity = NA_real_, Specificity = NA_real_, MCC = NA_real_, Brier = NA_real_))
  pred <- as.integer(score >= threshold)
  tp <- sum(pred == 1L & y_true == 1L); tn <- sum(pred == 0L & y_true == 0L)
  fp <- sum(pred == 1L & y_true == 0L); fn <- sum(pred == 0L & y_true == 1L)
  sens <- if (tp + fn > 0L) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0L) tn / (tn + fp) else NA_real_
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (is.finite(den) && den > 0) (tp * tn - fp * fn) / den else NA_real_
  ba <- if (is.finite(sens) && is.finite(spec)) mean(c(sens, spec)) else NA_real_
  brier <- if (is_probability) mean((score - y_true)^2) else NA_real_
  c(AUC = auc_rank(y_true, score), BalancedAccuracy = ba,
    Sensitivity = sens, Specificity = spec, MCC = mcc, Brier = brier)
}

# Return fold-specific RNA logCPM. Training TMM factors are estimated using ONLY
# training donors. Each validation/test donor is normalized against one fixed
# reference donor selected from the training partition, then placed on the same
# training normalization scale.
fold_rna_logcpm <- function(train_counts, test_counts) {
  train_counts <- as.matrix(train_counts); test_counts <- as.matrix(test_counts)
  suppressWarnings(storage.mode(train_counts) <- "numeric")
  suppressWarnings(storage.mode(test_counts) <- "numeric")
  train_counts[!is.finite(train_counts) | is.na(train_counts)] <- 0
  test_counts[!is.finite(test_counts) | is.na(test_counts)] <- 0
  train_counts <- round(pmax(train_counts, 0))
  test_counts <- round(pmax(test_counts, 0))
  keep_gene <- rowSums(train_counts) > 0
  train_counts <- train_counts[keep_gene, , drop = FALSE]
  test_counts <- test_counts[keep_gene, , drop = FALSE]
  if (nrow(train_counts) < 100L) stop("Too few nonzero RNA genes in inner/outer training split.")

  lib_train <- colSums(train_counts)
  lib_test <- if (ncol(test_counts)) colSums(test_counts) else numeric(0)
  if (any(!is.finite(lib_train) | lib_train <= 0)) stop("Zero RNA library size in training split.")
  if (length(lib_test) && any(!is.finite(lib_test) | lib_test <= 0)) {
    stop("Zero RNA library size in a held-out donor; cannot normalize that sample safely.")
  }
  ref_col <- which.min(abs(log(lib_train) - median(log(lib_train))))

  dge_train <- edgeR::DGEList(counts = train_counts)
  dge_train <- edgeR::calcNormFactors(dge_train, method = "TMM", refColumn = ref_col)
  nf_train <- dge_train$samples$norm.factors
  if (any(!is.finite(nf_train) | nf_train <= 0)) stop("Invalid training TMM factors.")

  nf_test <- rep(NA_real_, ncol(test_counts))
  if (ncol(test_counts) > 0L) {
    ref_counts <- train_counts[, ref_col]
    for (j in seq_len(ncol(test_counts))) {
      pair_counts <- cbind(reference = ref_counts, heldout = test_counts[, j])
      pair_dge <- edgeR::DGEList(counts = pair_counts)
      pair_dge <- tryCatch(edgeR::calcNormFactors(pair_dge, method = "TMM", refColumn = 1L),
                           error = function(e) NULL)
      if (is.null(pair_dge)) {
        nf_test[j] <- nf_train[ref_col]
      } else {
        pair_nf <- pair_dge$samples$norm.factors
        relative <- pair_nf[2L] / pair_nf[1L]
        nf_test[j] <- nf_train[ref_col] * relative
      }
    }
  }
  nf_test[!is.finite(nf_test) | nf_test <= 0] <- nf_train[ref_col]

  effective_train <- lib_train * nf_train
  effective_test <- lib_test * nf_test
  cpm_train <- sweep(train_counts, 2L, effective_train, "/") * 1e6
  cpm_test <- if (ncol(test_counts)) sweep(test_counts, 2L, effective_test, "/") * 1e6 else
    matrix(numeric(0), nrow = nrow(train_counts), ncol = 0L,
           dimnames = list(rownames(train_counts), NULL))
  # log2(CPM + 1) is computed sample-wise with the training-derived TMM factors;
  # unlike a cohort-wide logCPM call, the prior is not recalculated using test libraries.
  list(train = log2(cpm_train + 1),
       test = if (ncol(test_counts)) log2(cpm_test + 1) else cpm_test,
       train_effective_library = effective_train,
       reference_training_donor = colnames(train_counts)[ref_col])
}

# Fold-specific preprocessing. Each layer is capped so combined models carry at
# most about 30 variables total; a single-layer model may use up to 20.
prep_block <- function(M, train_idx, test_idx, layer, max_features) {
  M <- as.matrix(M)
  suppressWarnings(storage.mode(M) <- "numeric")
  if (is.null(rownames(M))) rownames(M) <- paste0(layer, "_feature_", seq_len(nrow(M)))
  tr_raw <- M[, train_idx, drop = FALSE]
  te_raw <- M[, test_idx, drop = FALSE]
  if (any(is.infinite(tr_raw)) || any(is.infinite(te_raw))) stop("Infinite values in ", layer)

  if (layer == "RNA") {
    norm <- fold_rna_logcpm(tr_raw, te_raw)
    tr <- norm$train
    te <- norm$test
    # Expression filter learned from training donors only: CPM >= 1 in at least
    # 20% of training donors, with a minimum of two donors.
    effective_lib <- norm$train_effective_library
    raw_cpm_train <- sweep(tr_raw[rownames(tr), , drop = FALSE], 2L,
                           effective_lib, "/") * 1e6
    min_detect <- max(2L, ceiling(0.20 * ncol(tr)))
    expressed <- rowSums(raw_cpm_train >= 1, na.rm = TRUE) >= min_detect
    tr <- tr[expressed, , drop = FALSE]
    te <- te[expressed, , drop = FALSE]
  } else {
    tr <- pmax(tr_raw, 0)
    te <- pmax(te_raw, 0)
    # Peak-gene features must be detected in at least 20% of training donors.
    min_detect <- max(2L, ceiling(0.20 * ncol(tr)))
    detected <- rowSums(tr > 0, na.rm = TRUE) >= min_detect
    tr <- log1p(tr[detected, , drop = FALSE])
    te <- log1p(te[detected, , drop = FALSE])
  }

  if (nrow(tr) == 0L) stop("No features remain after training-only prevalence filter for ", layer)

  # Training-median imputation only; test samples receive training medians.
  med <- matrixStats::rowMedians(tr, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  bad_tr <- which(!is.finite(tr), arr.ind = TRUE)
  if (nrow(bad_tr)) tr[bad_tr] <- med[bad_tr[, 1L]]
  bad_te <- which(!is.finite(te), arr.ind = TRUE)
  if (nrow(bad_te)) te[bad_te] <- med[bad_te[, 1L]]

  vars <- matrixStats::rowVars(tr, na.rm = TRUE)
  vars[!is.finite(vars)] <- 0
  eligible <- which(vars > 1e-10)
  if (!length(eligible)) stop("No variable features remain for ", layer)
  ord <- eligible[order(vars[eligible], decreasing = TRUE)]
  ord <- head(ord, max_features)
  tr <- tr[ord, , drop = FALSE]
  te <- te[ord, , drop = FALSE]

  mu <- rowMeans(tr)
  sig <- matrixStats::rowSds(tr)
  keep <- is.finite(sig) & sig > 1e-8
  tr <- tr[keep, , drop = FALSE]
  te <- te[keep, , drop = FALSE]
  mu <- mu[keep]; sig <- sig[keep]
  if (nrow(tr) == 0L) stop("All selected features were constant in training data for ", layer)

  tr <- sweep(sweep(tr, 1L, mu, "-"), 1L, sig, "/")
  te <- sweep(sweep(te, 1L, mu, "-"), 1L, sig, "/")
  colnames(tr) <- colnames(M)[train_idx]
  colnames(te) <- colnames(M)[test_idx]
  list(x_train = t(tr), x_test = t(te),
       features = paste0(layer, "::", rownames(tr)),
       n_after_filter = nrow(tr))
}

# Stratified inner and outer folds are the only units used for learning.
# The cache key includes the exact train/test donor indices and feature cap.
layer_cache <- new.env(parent = emptyenv())
get_prepped_block <- function(block, train_idx, test_idx, cap) {
  key <- paste(block, paste(train_idx, collapse = ","), paste(test_idx, collapse = ","), cap, sep = "__")
  if (exists(key, envir = layer_cache, inherits = FALSE)) return(get(key, envir = layer_cache, inherits = FALSE))
  value <- prep_block(obj[[block]], train_idx, test_idx, block, cap)
  assign(key, value, envir = layer_cache)
  value
}

make_design <- function(block_names, train_idx, test_idx) {
  cap <- min(max_features_per_layer, max(5L, floor(max_total_features / length(block_names))))
  parts <- lapply(block_names, function(b) get_prepped_block(b, train_idx, test_idx, cap))
  x_train <- do.call(cbind, lapply(parts, `[[`, "x_train"))
  x_test <- do.call(cbind, lapply(parts, `[[`, "x_test"))
  feature_names <- make.unique(unlist(lapply(parts, `[[`, "features"), use.names = FALSE))
  if (is.null(dim(x_train)) || ncol(x_train) == 0L ||
      ncol(x_train) != ncol(x_test) || length(feature_names) != ncol(x_train)) {
    stop("Empty or inconsistent design matrix for model: ", paste(block_names, collapse = "+"))
  }
  colnames(x_train) <- feature_names; colnames(x_test) <- feature_names
  list(x_train = x_train, x_test = x_test, features = feature_names,
       feature_counts = data.frame(modality = block_names,
         n_features = vapply(parts, function(z) length(z$features), integer(1))))
}

ratio_grid <- exp(seq(log(1), log(min_ratio), length.out = nlambda))
fit_en <- function(x_train, y_train, x_new, context) {
  captured <- character()
  fit <- withCallingHandlers(
    glmnet::glmnet(x = x_train, y = y_train, family = "binomial", alpha = alpha,
      nlambda = nlambda, lambda.min.ratio = min_ratio, standardize = FALSE,
      maxit = 100000),
    warning = function(w) {
      captured <<- c(captured, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(captured)) for (w in unique(captured)) record_warning(context, w)
  pred <- vapply(ratio_grid, function(r) {
    as.numeric(predict(fit, newx = x_new, s = max(fit$lambda) * r, type = "response"))
  }, numeric(nrow(x_new)))
  if (is.null(dim(pred))) pred <- matrix(pred, ncol = length(ratio_grid))
  list(fit = fit, pred_grid = pred, ratios = fit$lambda / max(fit$lambda))
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

all_pred <- list(); all_metrics <- list(); all_selected <- list()
feature_sets <- list(); tuning_rows <- list(); feature_qc_rows <- list()
completed_model_folds <- list(); warning_log <- data.frame(context = character(), warning = character())
warning_env <- new.env(parent = emptyenv())
record_warning <- function(context, message) {
  warning_env$rows <- rbind(warning_env$rows,
    data.frame(context = as.character(context), warning = as.character(message), stringsAsFactors = FALSE))
}
warning_env$rows <- warning_log

pred_counter <- 0L; selected_counter <- 0L; feature_set_counter <- 0L
completed_fold_counter <- 0L; tuning_counter <- 0L; feature_qc_counter <- 0L

for (rep_id in seq_len(repeats)) {
  outer_fold <- make_folds(y, outer_k, 20261117 + rep_id)
  rep_pred <- list()
  for (fold_id in seq_len(outer_k)) {
    rm(list = ls(envir = layer_cache, all.names = TRUE), envir = layer_cache)
    train_idx <- which(outer_fold != fold_id)
    test_idx <- which(outer_fold == fold_id)
    y_train <- y[train_idx]
    inner_fold <- make_folds(y_train, inner_k, 20261117 + rep_id * 100L + fold_id)

    for (model_name in names(sets)) {
      blocks <- sets[[model_name]]
      # Every row of inner_oof is predicted by a model that did not train on it.
      inner_oof <- matrix(NA_real_, nrow = length(train_idx), ncol = nlambda)
      inner_fold_auc <- matrix(NA_real_, nrow = inner_k, ncol = nlambda)
      inner_ok <- TRUE

      for (inner_id in seq_len(inner_k)) {
        inner_train_local <- which(inner_fold != inner_id)
        inner_valid_local <- which(inner_fold == inner_id)
        inner_train_idx <- train_idx[inner_train_local]
        inner_valid_idx <- train_idx[inner_valid_local]
        if (length(unique(y[inner_train_idx])) < 2L || length(unique(y[inner_valid_idx])) < 2L) {
          record_warning(paste("inner", rep_id, fold_id, model_name, inner_id),
                         "Inner split had one class; nested tuning cannot evaluate this split.")
          inner_ok <- FALSE
          break
        }
        design_inner <- tryCatch(make_design(blocks, inner_train_idx, inner_valid_idx),
          error = function(e) {
            record_warning(paste("inner", rep_id, fold_id, model_name, inner_id), conditionMessage(e)); NULL
          })
        if (is.null(design_inner)) { inner_ok <- FALSE; break }
        fit_inner <- tryCatch(fit_en(design_inner$x_train, y[inner_train_idx], design_inner$x_test,
                                     paste("inner", rep_id, fold_id, model_name, inner_id)),
          error = function(e) {
            record_warning(paste("inner", rep_id, fold_id, model_name, inner_id), conditionMessage(e)); NULL
          })
        if (is.null(fit_inner)) { inner_ok <- FALSE; break }
        inner_oof[inner_valid_local, ] <- fit_inner$pred_grid
        inner_fold_auc[inner_id, ] <- vapply(seq_len(nlambda), function(j) {
          auc_rank(y[inner_valid_idx], fit_inner$pred_grid[, j])
        }, numeric(1))
      }

      if (!inner_ok) {
        record_warning(paste("outer", rep_id, fold_id, model_name),
                       "Model/fold skipped because an inner split could not be completed.")
        next
      }
      pooled_inner_auc <- vapply(seq_len(nlambda), function(j) auc_rank(y_train, inner_oof[, j]), numeric(1))
      if (!any(is.finite(pooled_inner_auc))) {
        record_warning(paste("outer", rep_id, fold_id, model_name), "No valid pooled inner OOF AUC for lambda selection.")
        next
      }
      best_auc <- max(pooled_inner_auc, na.rm = TRUE)
      eligible <- which(is.finite(pooled_inner_auc) & pooled_inner_auc >= best_auc - inner_auc_tolerance)
      # Among near-optimal inner OOF scores, prefer the strongest regularization.
      selected_ratio <- max(ratio_grid[eligible])
      selected_idx <- which.min(abs(log(ratio_grid) - log(selected_ratio)))

      design_outer <- tryCatch(make_design(blocks, train_idx, test_idx),
        error = function(e) {
          record_warning(paste("outer", rep_id, fold_id, model_name), conditionMessage(e)); NULL
        })
      if (is.null(design_outer)) next
      fit_outer <- tryCatch(fit_en(design_outer$x_train, y_train, design_outer$x_test,
                                   paste("outer", rep_id, fold_id, model_name)),
        error = function(e) {
          record_warning(paste("outer", rep_id, fold_id, model_name), conditionMessage(e)); NULL
        })
      if (is.null(fit_outer)) next
      lambda_target <- max(fit_outer$fit$lambda) * selected_ratio
      lambda_idx <- which.min(abs(log(fit_outer$fit$lambda) - log(lambda_target)))
      lambda <- fit_outer$fit$lambda[lambda_idx]
      probability <- as.numeric(predict(fit_outer$fit, newx = design_outer$x_test,
                                        s = lambda, type = "response"))
      if (length(probability) != length(test_idx) || any(!is.finite(probability))) {
        record_warning(paste("outer", rep_id, fold_id, model_name), "Invalid outer-fold probabilities; fold skipped.")
        next
      }
      metric <- calc_metrics(y[test_idx], probability)
      all_metrics[[length(all_metrics) + 1L]] <- data.frame(
        repetition = rep_id, fold = fold_id, model = model_name,
        AUC = unname(metric["AUC"]), BalancedAccuracy = unname(metric["BalancedAccuracy"]),
        Sensitivity = unname(metric["Sensitivity"]), Specificity = unname(metric["Specificity"]),
        MCC = unname(metric["MCC"]), Brier = unname(metric["Brier"]),
        n_test = length(test_idx), stringsAsFactors = FALSE
      )
      pred_counter <- pred_counter + 1L
      pred_row <- data.frame(
        repetition = rep_id, fold = fold_id, model = model_name,
        donor = manifest$donor[test_idx], group = as.character(manifest$group[test_idx]),
        truth = y[test_idx], probability = probability, decision_score = NA_real_,
        score_type = "probability", threshold = 0.5, stringsAsFactors = FALSE
      )
      rep_pred[[length(rep_pred) + 1L]] <- pred_row

      co <- as.matrix(stats::coef(fit_outer$fit, s = lambda))[, 1L]
      nz <- which(is.finite(co) & abs(co) > 1e-10 & names(co) != "(Intercept)")
      selected_names <- if (length(nz)) names(co)[nz] else character()
      feature_set_counter <- feature_set_counter + 1L
      feature_sets[[feature_set_counter]] <- list(model = model_name, repetition = rep_id,
                                                  fold = fold_id, feature_id = selected_names)
      if (length(nz)) {
        selected_counter <- selected_counter + 1L
        all_selected[[selected_counter]] <- data.frame(
          model = model_name, feature_id = names(co)[nz], coefficient = as.numeric(co[nz]),
          repetition = rep_id, fold = fold_id, stringsAsFactors = FALSE
        )
      }
      completed_fold_counter <- completed_fold_counter + 1L
      completed_model_folds[[completed_fold_counter]] <- data.frame(
        model = model_name, repetition = rep_id, fold = fold_id, stringsAsFactors = FALSE
      )
      tuning_counter <- tuning_counter + 1L
      tuning_rows[[tuning_counter]] <- data.frame(
        repetition = rep_id, outer_fold = fold_id, model = model_name,
        alpha = alpha, selected_lambda = lambda, selected_lambda_ratio = selected_ratio,
        best_inner_pooled_oof_auc = best_auc,
        selected_inner_pooled_oof_auc = pooled_inner_auc[selected_idx],
        inner_oof_donors = sum(complete.cases(inner_oof[, selected_idx, drop = FALSE])),
        stringsAsFactors = FALSE
      )
      feature_qc_counter <- feature_qc_counter + 1L
      feature_qc_rows[[feature_qc_counter]] <- data.frame(
        repetition = rep_id, outer_fold = fold_id, model = model_name,
        n_features_total = ncol(design_outer$x_train),
        feature_counts = paste(paste(design_outer$feature_counts$modality,
          design_outer$feature_counts$n_features, sep = "="), collapse = ";"),
        n_nonzero_coefficients = length(nz), stringsAsFactors = FALSE
      )
      cat("RA rep", rep_id, "fold", fold_id, model_name,
          "inner AUC", round(best_auc, 3), "outer AUC", round(metric["AUC"], 3),
          "features", ncol(design_outer$x_train), "selected", length(nz), "\n")
    }
  }

  if (!length(rep_pred)) stop("No RA models produced predictions in repetition ", rep_id)
  rp <- do.call(rbind, rep_pred)
  pred_rep_metrics <- do.call(rbind, lapply(split(rp, rp$model), function(z) {
    mm <- calc_metrics(z$truth, z$probability)
    data.frame(repetition = rep_id, fold = "pooled_OOF", model = z$model[1],
      AUC = unname(mm["AUC"]), BalancedAccuracy = unname(mm["BalancedAccuracy"]),
      Sensitivity = unname(mm["Sensitivity"]), Specificity = unname(mm["Specificity"]),
      MCC = unname(mm["MCC"]), Brier = unname(mm["Brier"]), n_test = nrow(z),
      stringsAsFactors = FALSE)
  }))
  all_metrics <- c(all_metrics, lapply(seq_len(nrow(pred_rep_metrics)), function(i) pred_rep_metrics[i, , drop = FALSE]))
  all_pred <- c(all_pred, list(rp))
  utils::write.csv(do.call(rbind, all_pred), file.path(out_dir, "ra_oof_predictions_partial.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, all_metrics), file.path(out_dir, "ra_metrics_partial.csv"), row.names = FALSE)
  cat("Checkpoint saved after RA repetition", rep_id, "of", repeats, "\n")
}

# Fixed linear SVM sensitivity analysis. C=0.5 is pre-specified; no tuning on
# outer-test donors. Decision scores are stored as scores, NOT probabilities.
if (run_svm) {
  svm_specs <- list(RNA_SVM = "RNA", ALL_SVM = required_blocks)
  for (rep_id in seq_len(repeats)) {
    outer_fold <- make_folds(y, outer_k, 20261117 + rep_id)
    svm_pred <- list()
    for (fold_id in seq_len(outer_k)) {
      rm(list = ls(envir = layer_cache, all.names = TRUE), envir = layer_cache)
      train_idx <- which(outer_fold != fold_id); test_idx <- which(outer_fold == fold_id)
      for (nm in names(svm_specs)) {
        des <- tryCatch(make_design(svm_specs[[nm]], train_idx, test_idx), error = function(e) {
          record_warning(paste("SVM", rep_id, fold_id, nm), conditionMessage(e)); NULL
        })
        if (is.null(des)) next
        fit <- tryCatch(e1071::svm(x = des$x_train, y = factor(y[train_idx], levels = c(0, 1)),
          kernel = "linear", cost = 0.5, scale = FALSE, probability = FALSE), error = function(e) {
            record_warning(paste("SVM", rep_id, fold_id, nm), conditionMessage(e)); NULL
          })
        if (is.null(fit)) next
        train_pred <- predict(fit, des$x_train, decision.values = TRUE)
        test_pred <- predict(fit, des$x_test, decision.values = TRUE)
        train_score <- as.numeric(attr(train_pred, "decision.values"))
        test_score <- as.numeric(attr(test_pred, "decision.values"))
        if (length(train_score) != length(train_idx) || length(test_score) != length(test_idx) ||
            any(!is.finite(train_score)) || any(!is.finite(test_score))) {
          record_warning(paste("SVM", rep_id, fold_id, nm), "Invalid SVM decision scores; fold skipped.")
          next
        }
        # SVM sign convention is resolved with training labels only.
        if (mean(train_score[y[train_idx] == 1L]) < mean(train_score[y[train_idx] == 0L])) {
          train_score <- -train_score; test_score <- -test_score
        }
        pred_obj <- data.frame(repetition = rep_id, fold = fold_id, model = nm,
          donor = manifest$donor[test_idx], group = as.character(manifest$group[test_idx]),
          truth = y[test_idx], probability = NA_real_, decision_score = test_score,
          score_type = "decision_score", threshold = 0, stringsAsFactors = FALSE)
        svm_pred[[length(svm_pred) + 1L]] <- pred_obj
        sm <- calc_metrics(y[test_idx], test_score, threshold = 0, is_probability = FALSE)
        all_metrics[[length(all_metrics) + 1L]] <- data.frame(
          repetition = rep_id, fold = fold_id, model = nm,
          AUC = unname(sm["AUC"]), BalancedAccuracy = unname(sm["BalancedAccuracy"]),
          Sensitivity = unname(sm["Sensitivity"]), Specificity = unname(sm["Specificity"]),
          MCC = unname(sm["MCC"]), Brier = NA_real_, n_test = length(test_idx),
          stringsAsFactors = FALSE)
      }
    }
    if (length(svm_pred)) {
      sp <- do.call(rbind, svm_pred)
      for (nm in unique(sp$model)) {
        z <- sp[sp$model == nm, , drop = FALSE]
        sm <- calc_metrics(z$truth, z$decision_score, threshold = 0, is_probability = FALSE)
        all_metrics[[length(all_metrics) + 1L]] <- data.frame(
          repetition = rep_id, fold = "pooled_OOF", model = nm,
          AUC = unname(sm["AUC"]), BalancedAccuracy = unname(sm["BalancedAccuracy"]),
          Sensitivity = unname(sm["Sensitivity"]), Specificity = unname(sm["Specificity"]),
          MCC = unname(sm["MCC"]), Brier = NA_real_, n_test = nrow(z),
          stringsAsFactors = FALSE)
      }
      all_pred[[length(all_pred) + 1L]] <- sp
    }
  }
} else if (isTRUE(a$ra_run_svm_sensitivity)) {
  record_warning("SVM", "e1071 is not installed; SVM sensitivity analysis was skipped.")
}

if (!length(all_pred)) stop("No RA predictions were generated.")
pred <- do.call(rbind, all_pred)
met <- do.call(rbind, all_metrics)
# Pooled outer-fold rows are used for headline metrics; fold-level values remain
# available for diagnostics. Repeated predictions are kept separate by repetition.
pooled <- met[met$fold == "pooled_OOF", , drop = FALSE]
metric_names <- c("AUC", "BalancedAccuracy", "Sensitivity", "Specificity", "MCC", "Brier")
metric_mean <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
metric_sd <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)

summary <- do.call(rbind, lapply(split(pooled, pooled$model), function(z) {
  row <- data.frame(model = z$model[1], n_repetitions = nrow(z), stringsAsFactors = FALSE)
  for (nm in metric_names) {
    safe_nm <- tolower(gsub("([a-z])([A-Z])", "\\1_\\2", nm))
    row[[paste0("mean_", safe_nm)]] <- metric_mean(z[[nm]])
    row[[paste0("sd_", safe_nm)]] <- metric_sd(z[[nm]])
  }
  row
}))

repeat_summary <- do.call(rbind, lapply(split(pooled, paste(pooled$repetition, pooled$model, sep = "__")), function(z) {
  row <- data.frame(repetition = z$repetition[1], model = z$model[1], stringsAsFactors = FALSE)
  for (nm in metric_names) {
    safe_nm <- tolower(gsub("([a-z])([A-Z])", "\\1_\\2", nm))
    row[[paste0("mean_", safe_nm)]] <- metric_mean(z[[nm]])
  }
  row
}))

# Incremental contribution vs RNA-only, paired within each outer repetition.
rna_base <- pooled[pooled$model == "RNA", c("repetition", "AUC", "BalancedAccuracy"), drop = FALSE]
if (nrow(rna_base)) {
  names(rna_base) <- c("repetition", "rna_auc", "rna_balanced_accuracy")
  contrib <- merge(pooled, rna_base, by = "repetition", all.x = TRUE)
  contrib$delta_auc_vs_RNA <- contrib$AUC - contrib$rna_auc
  contrib$delta_balanced_accuracy_vs_RNA <- contrib$BalancedAccuracy - contrib$rna_balanced_accuracy
  contrib_out <- do.call(rbind, lapply(split(contrib, contrib$model), function(z) data.frame(
    model = z$model[1], delta_auc_vs_RNA = metric_mean(z$delta_auc_vs_RNA),
    sd_delta_auc_vs_RNA = metric_sd(z$delta_auc_vs_RNA),
    delta_balanced_accuracy_vs_RNA = metric_mean(z$delta_balanced_accuracy_vs_RNA),
    sd_delta_balanced_accuracy_vs_RNA = metric_sd(z$delta_balanced_accuracy_vs_RNA),
    stringsAsFactors = FALSE)))
} else {
  contrib_out <- data.frame()
}

# Gene/feature selection frequency and pairwise Jaccard across successful outer folds.
selected_out <- data.frame(); stability_out <- data.frame(); jaccard_out <- data.frame()
if (length(all_selected)) {
  selected_out <- do.call(rbind, all_selected)
  selected_out <- selected_out[!duplicated(selected_out[c("model", "feature_id", "repetition", "fold")]), , drop = FALSE]
  selected_out$modality <- sub("::.*$", "", selected_out$feature_id)
  selected_out$gene_or_feature <- sub("^[^:]+::", "", selected_out$feature_id)
  selected_out$selected_count <- 1L
  stability_out <- aggregate(selected_count ~ model + feature_id + modality + gene_or_feature,
                             selected_out, sum)
  completed_df <- unique(do.call(rbind, completed_model_folds))
  denom <- aggregate(fold ~ model, completed_df, length)
  names(denom)[2L] <- "completed_outer_folds"
  stability_out <- merge(stability_out, denom, by = "model", all.x = TRUE)
  stability_out$selection_frequency <- stability_out$selected_count / stability_out$completed_outer_folds
  stability_out <- stability_out[order(stability_out$model, -stability_out$selection_frequency), ]
}
if (length(feature_sets)) {
  model_names <- unique(vapply(feature_sets, function(z) z$model, character(1)))
  for (nm in model_names) {
    items <- Filter(function(z) identical(z$model, nm), feature_sets)
    feature_lists <- lapply(items, function(z) unique(z$feature_id[!is.na(z$feature_id) & nzchar(z$feature_id)]))
    if (length(feature_lists) >= 2L) {
      pairs <- utils::combn(seq_along(feature_lists), 2L)
      vals <- apply(pairs, 2L, function(ix) {
        aa <- feature_lists[[ix[1L]]]; bb <- feature_lists[[ix[2L]]]
        uu <- union(aa, bb)
        if (!length(uu)) NA_real_ else length(intersect(aa, bb)) / length(uu)
      })
      jaccard_out <- rbind(jaccard_out, data.frame(model = nm,
        comparisons = sum(is.finite(vals)), mean_jaccard = metric_mean(vals),
        sd_jaccard = metric_sd(vals), stringsAsFactors = FALSE))
    }
  }
}

utils::write.csv(pred, file.path(out_dir, "ra_oof_predictions.csv"), row.names = FALSE)
utils::write.csv(met, file.path(out_dir, "ra_metrics_by_oof_repeat.csv"), row.names = FALSE)
utils::write.csv(repeat_summary, file.path(out_dir, "ra_metrics_by_repeat.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "ra_model_metrics.csv"), row.names = FALSE)
utils::write.csv(contrib_out, file.path(out_dir, "ra_layer_contribution.csv"), row.names = FALSE)
if (nrow(selected_out)) utils::write.csv(selected_out, file.path(out_dir, "ra_selected_features_by_fold.csv"), row.names = FALSE)
if (nrow(stability_out)) utils::write.csv(stability_out, file.path(out_dir, "ra_stable_features.csv"), row.names = FALSE)
if (nrow(jaccard_out)) utils::write.csv(jaccard_out, file.path(out_dir, "ra_stability_jaccard.csv"), row.names = FALSE)
if (length(tuning_rows)) utils::write.csv(do.call(rbind, tuning_rows), file.path(out_dir, "ra_inner_tuning.csv"), row.names = FALSE)
if (length(feature_qc_rows)) utils::write.csv(do.call(rbind, feature_qc_rows), file.path(out_dir, "ra_feature_qc_by_fold.csv"), row.names = FALSE)
if (nrow(warning_env$rows)) utils::write.csv(warning_env$rows, file.path(out_dir, "ra_model_warnings.csv"), row.names = FALSE)

qc_manifest <- data.frame(
  donors = length(common), RA = sum(y == 1L), OA = sum(y == 0L),
  outer_repeats = repeats, outer_folds = outer_k, inner_folds = inner_k,
  max_features_per_layer = max_features_per_layer, max_total_features_per_model = max_total_features,
  nlambda = nlambda, lambda_min_ratio = min_ratio, alpha = alpha,
  inner_auc_tolerance = inner_auc_tolerance, svm_sensitivity_run = run_svm,
  preprocessing = paste(
    "Raw RNA counts; TMM factors and logCPM estimated within each training split;",
    "held-out RNA libraries normalized relative to one fixed training reference;",
    "RNA expression prevalence filter, chromatin detection filter, imputation,",
    "variance ranking, and z-scoring fitted from training donors only."
  ),
  validation = paste(
    "Donor-stratified repeated nested CV; inner tuning uses pooled inner out-of-fold AUC;",
    "outer-test donors are excluded from preprocessing and model tuning."
  ),
  caveat = paste(
    "Only 20 donors are available (10 RA, 10 OA). Metrics are exploratory and can be unstable;",
    "repeated CV does not increase the number of independent donors. External validation is required."
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(qc_manifest, file.path(out_dir, "ra_qc_model_manifest.csv"), row.names = FALSE)

# Remove transient checkpoints only after final outputs were written successfully.
unlink(file.path(out_dir, c("ra_oof_predictions_partial.csv", "ra_metrics_partial.csv")))
cat("\nRA leakage-conscious benchmark completed. Donors:", length(common),
    "| predictions:", nrow(pred), "| metric rows:", nrow(met), "\n")
cat("Elastic Net configurations:", length(sets), "| SVM sensitivity:", run_svm,
    "| warnings:", nrow(warning_env$rows), "\n")
cat("Outputs written to:", out_dir, "\n")
