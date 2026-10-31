# Install only missing CRAN/Bioconductor dependencies for MOSAIC-3.
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=3600)
cran <- c("jsonlite","data.table","survival","glmnet","ranger","e1071","matrixStats")
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing,dependencies=TRUE)
if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
bioc <- c("edgeR","GenomicRanges","GenomicFeatures","rtracklayer","AnnotationDbi",
          "org.Hs.eg.db","GO.db","TxDb.Hsapiens.UCSC.hg19.knownGene")
missing_bioc <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc,ask=FALSE,update=FALSE)
cat("MOSAIC-3 dependencies verified.\n")
