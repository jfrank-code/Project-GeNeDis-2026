# MOSAIC-3 | Prepare matched RA/OA synovial fibroblast multi-omics from GEO.
# Uses processed RNA counts + per-sample ATAC BED + H3K27ac/H3K4me3 peaks.
# No WGBS, FASTQ/SRA, or full GEO TAR archive is downloaded or processed.
options(stringsAsFactors = FALSE, timeout = 3600)
for (p in c("data.table","edgeR","GenomicRanges","GenomicFeatures","rtracklayer",
            "AnnotationDbi","org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene")) {
  if (!requireNamespace(p, quietly=TRUE)) stop("Missing R package: ", p, ". Run Rscript src/scripts/R/00_install_dependencies.R")
}
root <- normalizePath(".", winslash="/", mustWork=TRUE)
raw <- file.path(root,"data","raw","RA")
proc <- file.path(root,"data","processed","RA")
dir.create(proc, recursive=TRUE, showWarnings=FALSE)
man <- read.csv(file.path(raw,"geo_download_manifest.csv"), stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("path","assay","donor") %in% names(man))) stop("Manifest GEO incompleto.")

# GEO's RNA matrix includes exactly 10 RA and 10 OA RNA donors. Resolve by donor
# labels in column names, never by assuming column order.
rna_file <- file.path(raw,"RNA","GSE112656_Counts_matix.txt.gz")
if (!file.exists(rna_file)) stop("Matrix not found: RNA GEO: ",rna_file)
rna_dt <- data.table::fread(rna_file, data.table=FALSE, check.names=FALSE)
if (ncol(rna_dt)<3) stop("GEO RNA matrix is not tabular.")
feature_id <- as.character(rna_dt[[1]])
rna_dt[[1]] <- NULL
rna_mat <- as.matrix(data.frame(lapply(rna_dt,function(z)suppressWarnings(as.numeric(z))),check.names=FALSE))
rownames(rna_mat) <- feature_id
if (anyNA(rna_mat) || any(!is.finite(rna_mat)) || any(rna_mat<0)) stop("RNA counts contain invalid values.")
gsm_donor <- c(
  GSM3075549="RA_01",GSM3075550="RA_02",GSM3075551="RA_03",GSM3075552="RA_04",
  GSM3075553="RA_05",GSM3075554="RA_06",GSM3075555="RA_07",GSM3075556="RA_08",
  GSM3075557="RA_10",GSM3075558="RA_11",GSM3075559="OA_02",GSM3075560="OA_03",
  GSM3075561="OA_04",GSM3075562="OA_05",GSM3075563="OA_06",GSM3075564="OA_07",
  GSM3075565="OA_08",GSM3075566="OA_09",GSM3075567="OA_10",GSM3075568="OA_11")
extract_donor <- function(x) {
  ux <- toupper(as.character(x))
  m <- regmatches(ux,regexpr("(RA|OA)_[0-9]{2}",ux))
  if (length(m)>0 && nzchar(m)) return(m)
  gsm <- regmatches(ux,regexpr("GSM[0-9]+",ux))
  if (length(gsm)>0 && nzchar(gsm) && gsm %in% names(gsm_donor)) return(unname(gsm_donor[gsm]))
  NA_character_
}
col_donors <- vapply(colnames(rna_mat),extract_donor,character(1))
if (anyNA(col_donors) || anyDuplicated(col_donors)) stop("Could not resolve RNA columns to unique donor IDs. Columns: ",paste(head(colnames(rna_mat),10),collapse=", "))
colnames(rna_mat) <- col_donors
# Map feature IDs to HGNC symbols; unmapped IDs are retained with a stable prefix.
map_symbols <- function(ids) {
  x <- sub("\\.[0-9]+$","",as.character(ids)); out <- x
  looks_ens <- grepl("^ENSG[0-9]+$",x)
  looks_entrez <- grepl("^[0-9]+$",x)
  if (any(looks_ens)) {
    m <- tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=unique(x[looks_ens]),
        keytype="ENSEMBL",column="SYMBOL",multiVals="first"),error=function(e)NULL)
    if (!is.null(m)) { z <- unname(m[match(x,names(m))]); ok <- !is.na(z)&nzchar(z); out[ok] <- z[ok] }
  }
  if (any(looks_entrez)) {
    m <- tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=unique(x[looks_entrez]),
        keytype="ENTREZID",column="SYMBOL",multiVals="first"),error=function(e)NULL)
    if (!is.null(m)) { z <- unname(m[match(x,names(m))]); ok <- !is.na(z)&nzchar(z); out[ok] <- z[ok] }
  }
  # Already-symbol identifiers are retained; unmapped Ensembl/Entrez IDs remain traceable.
  out[!nzchar(out)|is.na(out)] <- paste0("UNMAPPED_",which(!nzchar(out)|is.na(out)))
  out
}

rna_symbols <- map_symbols(rownames(rna_mat))
# Collapse multiple transcript/Ensembl rows mapping to the same HGNC symbol.
rna_mat <- rowsum(rna_mat, group=rna_symbols, reorder=FALSE)
# Convert raw counts to logCPM. The same transform is used only within the RA block;
# ML filtering/imputation/scaling will still be estimated in each training fold.
y <- edgeR::DGEList(counts=rna_mat)
y <- edgeR::calcNormFactors(y, method="TMM")
rna_logcpm <- edgeR::cpm(y, log=TRUE, prior.count=1)

# hg19 promoters: upstream 2 kb, downstream 500 bp relative to annotated TSS.
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
gene_gr <- GenomicFeatures::genes(txdb)
prom <- GenomicRanges::promoters(gene_gr, upstream=2000, downstream=500)
entrez_ids <- names(prom)
keep_prom <- !is.na(entrez_ids) & nzchar(entrez_ids)
prom <- prom[keep_prom]; entrez_ids <- entrez_ids[keep_prom]
prom_syms <- map_symbols(entrez_ids)

peak_to_gene <- function(path) {
  gr <- rtracklayer::import(path)
  if (length(gr)==0) return(setNames(numeric(0),character(0)))
  cols <- names(S4Vectors::mcols(gr))
  val <- if ("signalValue" %in% cols) as.numeric(S4Vectors::mcols(gr)$signalValue) else if ("score" %in% cols) as.numeric(S4Vectors::mcols(gr)$score) else rep(1,length(gr))
  val[!is.finite(val)|is.na(val)] <- 0
  hits <- GenomicRanges::findOverlaps(gr,prom,ignore.strand=TRUE)
  if (!length(hits)) return(setNames(numeric(0),character(0)))
  gene <- prom_syms[S4Vectors::subjectHits(hits)]
  agg <- tapply(val[S4Vectors::queryHits(hits)],gene,sum,na.rm=TRUE)
  setNames(as.numeric(agg),names(agg))
}
make_peak_block <- function(label) {
  z <- man[man$assay==label & man$donor %in% colnames(rna_logcpm),,drop=FALSE]
  if (!nrow(z)) stop("No files found for ",label," among RNA donors.")
  vals <- list()
  for (i in seq_len(nrow(z))) {
    path <- z$path[i]
    if (!file.exists(path)) stop("Processed peak file missing: ",path)
    vals[[z$donor[i]]] <- peak_to_gene(path)
  }
  genes <- sort(unique(unlist(lapply(vals,names),use.names=FALSE)))
  M <- matrix(0,nrow=length(genes),ncol=length(vals),dimnames=list(genes,names(vals)))
  for (d in names(vals)) if (length(vals[[d]])) M[names(vals[[d]]),d] <- vals[[d]]
  # Collapse duplicate mapped genes and drop non-varying features only later within folds.
  M
}

atac <- make_peak_block("ATAC")
h27 <- make_peak_block("H3K27ac")
h4 <- make_peak_block("H3K4me3")
common <- Reduce(intersect,list(colnames(rna_logcpm),colnames(atac),colnames(h27),colnames(h4)))
if (length(common)<12) stop("Insufficient complete RNA+ATAC+H3K27ac+H3K4me3 intersection: ",length(common)," donors. Check the manifest.")
# Stable, explicit order; all blocks now represent the exact same biological donors.
common <- sort(common)
modalities <- list(RNA=rna_logcpm[,common,drop=FALSE], ATAC=atac[,common,drop=FALSE],
                   H3K27ac=h27[,common,drop=FALSE], H3K4me3=h4[,common,drop=FALSE])
manifest <- data.frame(donor=common,group=ifelse(grepl("^RA_",common),"RA","OA"),
                       RNA=common %in% colnames(rna_logcpm),ATAC=common %in% colnames(atac),
                       H3K27ac=common %in% colnames(h27),H3K4me3=common %in% colnames(h4),
                       stringsAsFactors=FALSE)
availability <- data.frame(modality=names(modalities),features=vapply(modalities,nrow,integer(1)),
                           donors_available=c(ncol(rna_logcpm),ncol(atac),ncol(h27),ncol(h4)),
                           donors_complete=length(common),RA=sum(manifest$group=="RA"),OA=sum(manifest$group=="OA"))
saveRDS(modalities,file.path(proc,"ra_modalities.rds"),compress="xz")
write.csv(manifest,file.path(proc,"ra_manifest.csv"),row.names=FALSE)
write.csv(availability,file.path(proc,"ra_layer_qc.csv"),row.names=FALSE)
write.csv(man,file.path(proc,"geo_download_manifest.csv"),row.names=FALSE)
cat("RA prepared with",length(common),"matched donors:",sum(manifest$group=="RA"),"RA y",sum(manifest$group=="OA"),"OA.\n")
cat("Features RNA/ATAC/H3K27ac/H3K4me3:",paste(availability$features,collapse=" / "),"\n")
cat("WGBS excluded. Individual GEO files are retained for traceability.\n")
