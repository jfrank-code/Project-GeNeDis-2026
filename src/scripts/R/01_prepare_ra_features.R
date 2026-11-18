# MOSAIC-3 | Leakage-conscious preparation of paired RA/OA synovial multi-omics
#
# IMPORTANT DESIGN CHANGES
# 1) RNA raw counts are saved; TMM/logCPM is NOT calculated here. It is estimated
#    separately inside every inner/outer training split by 03_run_ra_benchmark.
# 2) ATAC/H3K27ac/H3K4me3 use a fixed hg19 promoter-gene universe from TxDb,
#    rather than a feature universe defined by peaks observed across all donors.
# 3) Peak signal is normalized within each donor (per-million mapped promoter
#    signal). This uses only that donor's own assay values, not labels or other donors.
# 4) Donor matching, duplicate checks, and sample availability are explicit.
#
# This script prepares features only. It does not use RA/OA labels for feature
# selection and does not fit a predictive model.

options(stringsAsFactors = FALSE, timeout = 3600, warn = 1)

required <- c("data.table", "edgeR", "GenomicRanges", "GenomicFeatures",
              "rtracklayer", "AnnotationDbi", "org.Hs.eg.db",
              "TxDb.Hsapiens.UCSC.hg19.knownGene")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       ". Run Rscript src/scripts/R/00_install_dependencies.R")
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
raw <- file.path(root, "data", "raw", "RA")
proc <- file.path(root, "data", "processed", "RA")
dir.create(proc, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(raw, "geo_download_manifest.csv")
if (!file.exists(manifest_path)) stop("Missing GEO download manifest: ", manifest_path)
man <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("path", "assay", "donor") %in% names(man))) {
  stop("geo_download_manifest.csv must contain path, assay and donor columns.")
}
man$path <- as.character(man$path)
man$assay <- as.character(man$assay)
man$donor <- as.character(man$donor)

# GEO's processed RNA count matrix contains 10 RA and 10 OA donors.
rna_file <- file.path(raw, "RNA", "GSE112656_Counts_matix.txt.gz")
if (!file.exists(rna_file)) stop("RNA count matrix not found: ", rna_file)
rna_dt <- data.table::fread(rna_file, data.table = FALSE, check.names = FALSE)
if (ncol(rna_dt) < 3L) stop("GEO RNA matrix is not tabular or has fewer than two samples.")
feature_id <- as.character(rna_dt[[1L]])
rna_dt[[1L]] <- NULL
rna_mat <- as.matrix(data.frame(
  lapply(rna_dt, function(z) suppressWarnings(as.numeric(z))),
  check.names = FALSE
))
rownames(rna_mat) <- feature_id
if (anyNA(rna_mat) || any(!is.finite(rna_mat)) || any(rna_mat < 0)) {
  stop("RNA count matrix contains missing, non-finite or negative values.")
}
# Counts should be integer-like. Round only floating-point representation noise.
if (any(abs(rna_mat - round(rna_mat)) > 1e-6)) {
  stop("RNA matrix does not appear to contain raw integer counts. Check GEO source.")
}
rna_mat <- round(rna_mat)

gsm_donor <- c(
  GSM3075549 = "RA_01", GSM3075550 = "RA_02", GSM3075551 = "RA_03",
  GSM3075552 = "RA_04", GSM3075553 = "RA_05", GSM3075554 = "RA_06",
  GSM3075555 = "RA_07", GSM3075556 = "RA_08", GSM3075557 = "RA_10",
  GSM3075558 = "RA_11", GSM3075559 = "OA_02", GSM3075560 = "OA_03",
  GSM3075561 = "OA_04", GSM3075562 = "OA_05", GSM3075563 = "OA_06",
  GSM3075564 = "OA_07", GSM3075565 = "OA_08", GSM3075566 = "OA_09",
  GSM3075567 = "OA_10", GSM3075568 = "OA_11"
)

extract_donor <- function(x) {
  ux <- toupper(as.character(x))
  m <- regmatches(ux, regexpr("(RA|OA)_[0-9]{2}", ux))
  if (length(m) == 1L && !is.na(m) && nzchar(m)) return(m)
  gsm <- regmatches(ux, regexpr("GSM[0-9]+", ux))
  if (length(gsm) == 1L && !is.na(gsm) && nzchar(gsm) && gsm %in% names(gsm_donor)) {
    return(unname(gsm_donor[gsm]))
  }
  NA_character_
}

rna_donors <- vapply(colnames(rna_mat), extract_donor, character(1))
if (anyNA(rna_donors) || anyDuplicated(rna_donors)) {
  stop("Could not resolve RNA columns to unique donor IDs. Found: ",
       paste(colnames(rna_mat), collapse = ", "))
}
colnames(rna_mat) <- rna_donors

# Map Ensembl/Entrez identifiers to HGNC symbols and collapse duplicated symbols.
map_symbols <- function(ids) {
  x <- sub("\\.[0-9]+$", "", as.character(ids))
  out <- x
  is_ens <- grepl("^ENSG[0-9]+$", x)
  is_entrez <- grepl("^[0-9]+$", x)
  if (any(is_ens)) {
    m <- tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(x[is_ens]), keytype = "ENSEMBL", column = "SYMBOL",
        multiVals = "first"), error = function(e) NULL
    )
    if (!is.null(m)) {
      z <- unname(m[match(x, names(m))]); ok <- !is.na(z) & nzchar(z)
      out[ok] <- z[ok]
    }
  }
  if (any(is_entrez)) {
    m <- tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(x[is_entrez]), keytype = "ENTREZID", column = "SYMBOL",
        multiVals = "first"), error = function(e) NULL
    )
    if (!is.null(m)) {
      z <- unname(m[match(x, names(m))]); ok <- !is.na(z) & nzchar(z)
      out[ok] <- z[ok]
    }
  }
  bad <- is.na(out) | !nzchar(out)
  out[bad] <- paste0("UNMAPPED_", which(bad))
  out
}

rna_symbols <- map_symbols(rownames(rna_mat))
rna_counts <- rowsum(rna_mat, group = rna_symbols, reorder = FALSE)
if (anyDuplicated(rownames(rna_counts))) stop("RNA gene-symbol collapse left duplicated row names.")
if (any(colSums(rna_counts) <= 0)) stop("At least one RNA donor has zero total library size.")

# Fixed promoter definition: hg19, 2 kb upstream and 500 bp downstream of TSS.
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
gene_gr <- GenomicFeatures::genes(txdb)
prom <- GenomicRanges::promoters(gene_gr, upstream = 2000, downstream = 500)
entrez_ids <- names(prom)
valid_prom <- !is.na(entrez_ids) & nzchar(entrez_ids)
prom <- prom[valid_prom]
entrez_ids <- entrez_ids[valid_prom]
prom_symbols <- map_symbols(entrez_ids)
valid_symbols <- !is.na(prom_symbols) & nzchar(prom_symbols) & !grepl("^UNMAPPED_", prom_symbols)
prom <- prom[valid_symbols]
prom_symbols <- prom_symbols[valid_symbols]
# This universe is annotation-defined and therefore does not depend on which
# donors happen to contain peaks.
promoter_gene_universe <- sort(unique(prom_symbols))
if (length(promoter_gene_universe) < 1000L) {
  stop("Unexpectedly small mapped promoter universe: ", length(promoter_gene_universe))
}

peak_to_gene <- function(path) {
  gr <- rtracklayer::import(path)
  if (length(gr) == 0L) return(setNames(numeric(0), character(0)))
  mc <- S4Vectors::mcols(gr)
  cols <- names(mc)
  if ("signalValue" %in% cols) {
    val <- suppressWarnings(as.numeric(mc$signalValue))
  } else if ("score" %in% cols) {
    val <- suppressWarnings(as.numeric(mc$score))
  } else {
    val <- rep(1, length(gr))
  }
  val[!is.finite(val) | is.na(val)] <- 0
  val <- pmax(val, 0)
  hits <- GenomicRanges::findOverlaps(gr, prom, ignore.strand = TRUE)
  if (length(hits) == 0L) return(setNames(numeric(0), character(0)))
  gene <- prom_symbols[S4Vectors::subjectHits(hits)]
  qval <- val[S4Vectors::queryHits(hits)]
  ok <- !is.na(gene) & nzchar(gene) & is.finite(qval)
  if (!any(ok)) return(setNames(numeric(0), character(0)))
  agg <- tapply(qval[ok], gene[ok], sum, na.rm = TRUE)
  setNames(as.numeric(agg), names(agg))
}

resolve_manifest_path <- function(p) {
  p <- as.character(p)
  if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  candidate <- file.path(root, p)
  if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/"))
  stop("Peak file listed in GEO manifest does not exist: ", p)
}

make_peak_block <- function(label) {
  z <- man[man$assay == label & man$donor %in% colnames(rna_counts), , drop = FALSE]
  if (!nrow(z)) stop("No files found for assay ", label, " among RNA donors.")
  if (anyDuplicated(z$donor)) {
    dup <- unique(z$donor[duplicated(z$donor)])
    stop("More than one manifest file per donor for ", label, ": ", paste(dup, collapse = ", "))
  }
  missing_donors <- setdiff(colnames(rna_counts), z$donor)
  if (length(missing_donors)) {
    warning(label, " missing donors: ", paste(missing_donors, collapse = ", "))
  }
  vals <- vector("list", nrow(z))
  names(vals) <- z$donor
  for (i in seq_len(nrow(z))) {
    path <- resolve_manifest_path(z$path[i])
    vals[[z$donor[i]]] <- peak_to_gene(path)
  }
  donors <- names(vals)
  M <- matrix(0, nrow = length(promoter_gene_universe), ncol = length(donors),
              dimnames = list(promoter_gene_universe, donors))
  for (donor in donors) {
    v <- vals[[donor]]
    if (length(v)) {
      keep <- intersect(names(v), rownames(M))
      if (length(keep)) M[keep, donor] <- v[keep]
    }
  }
  # Per-donor total-signal normalization. This is sample-local and uses neither
  # diagnostic labels nor information from other donors.
  totals <- colSums(M)
  for (j in seq_len(ncol(M))) {
    if (is.finite(totals[j]) && totals[j] > 0) {
      M[, j] <- M[, j] / totals[j] * 1e6
    }
  }
  M
}

atac <- make_peak_block("ATAC")
h27 <- make_peak_block("H3K27ac")
h4 <- make_peak_block("H3K4me3")

common <- Reduce(intersect, list(colnames(rna_counts), colnames(atac),
                                 colnames(h27), colnames(h4)))
if (length(common) < 12L) {
  stop("Insufficient complete RNA+ATAC+H3K27ac+H3K4me3 intersection: ",
       length(common), " donors. Check the GEO manifest.")
}
common <- sort(common)
modalities <- list(
  RNA = rna_counts[, common, drop = FALSE],
  ATAC = atac[, common, drop = FALSE],
  H3K27ac = h27[, common, drop = FALSE],
  H3K4me3 = h4[, common, drop = FALSE]
)
manifest <- data.frame(
  donor = common,
  group = ifelse(grepl("^RA_", common), "RA", "OA"),
  RNA = common %in% colnames(rna_counts),
  ATAC = common %in% colnames(atac),
  H3K27ac = common %in% colnames(h27),
  H3K4me3 = common %in% colnames(h4),
  stringsAsFactors = FALSE
)
if (anyDuplicated(manifest$donor) || anyNA(manifest$group)) stop("RA donor manifest is invalid.")
if (sum(manifest$group == "RA") < 2L || sum(manifest$group == "OA") < 2L) {
  stop("Both RA and OA require at least two complete donors.")
}

availability <- data.frame(
  modality = names(modalities),
  features = vapply(modalities, nrow, integer(1)),
  donors_available = c(ncol(rna_counts), ncol(atac), ncol(h27), ncol(h4)),
  donors_complete = length(common),
  RA = sum(manifest$group == "RA"),
  OA = sum(manifest$group == "OA"),
  preprocessing = c("raw integer counts; TMM/logCPM estimated within CV training folds",
                    "fixed hg19 promoter universe; per-donor CPM-like signal",
                    "fixed hg19 promoter universe; per-donor CPM-like signal",
                    "fixed hg19 promoter universe; per-donor CPM-like signal"),
  stringsAsFactors = FALSE
)

saveRDS(modalities, file.path(proc, "ra_modalities.rds"), compress = "xz")
write.csv(manifest, file.path(proc, "ra_manifest.csv"), row.names = FALSE)
write.csv(availability, file.path(proc, "ra_layer_qc.csv"), row.names = FALSE)
# Preserve the acquisition manifest for traceability.
write.csv(man, file.path(proc, "geo_download_manifest.csv"), row.names = FALSE)

cat("RA prepared with", length(common), "matched donors:",
    sum(manifest$group == "RA"), "RA and", sum(manifest$group == "OA"), "OA.\n")
cat("Features RNA/ATAC/H3K27ac/H3K4me3:",
    paste(availability$features, collapse = " / "), "\n")
cat("RNA stored as raw counts; TMM/logCPM will be fitted inside CV training folds.\n")
cat("ATAC/ChIP promoter signals use a fixed hg19 gene universe and per-donor normalization.\n")
cat("No WGBS, FASTQ/SRA, or full GEO TAR archive was processed.\n")
