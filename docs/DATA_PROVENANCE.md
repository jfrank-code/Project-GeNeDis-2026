# Data provenance and source register

## NSCLC: TCGA-LUAD and TCGA-LUSC

**Primary repository:** UCSC Xena Public Data Hubs.  
**Host:** `https://tcga.xenahubs.net`  
**Survival hub:** `https://pancanatlas.xenahubs.net`

| Layer | LUAD object | LUSC object | Primary use |
|---|---|---|---|
| RNA-seq | `TCGA.LUAD.sampleMap/HiSeqV2` | `TCGA.LUSC.sampleMap/HiSeqV2` | Gene expression |
| CNV | `TCGA.LUAD.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes` | `TCGA.LUSC.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes` | Gene-level GISTIC2 copy number |
| Methylation | `TCGA.LUAD.sampleMap/HumanMethylation450` | `TCGA.LUSC.sampleMap/HumanMethylation450` | HumanMethylation450 beta values |
| Clinical | `TCGA.LUAD.sampleMap/LUAD_clinicalMatrix` | `TCGA.LUSC.sampleMap/LUSC_clinicalMatrix` | Age, sex, stage and case metadata |
| Survival | `survival/LUAD_survival.txt` | `survival/LUSC_survival.txt` | Overall-survival time and event |
| SNV (exploratory) | `mc3_gene_level/LUAD_mc3_gene_level.txt` | `mc3_gene_level/LUSC_mc3_gene_level.txt` | Gene-level MC3 mutation presence; excluded from primary model pending callability audit |

**Selection and harmonisation:** primary tumour sample type 01; one deterministic specimen per case (vial A preferred); intersection of RNA, CNV, HM450 and valid overall survival. Cohort sizes and endpoint counts are recorded in `results/NSCLC/nsclc_cohort_summary.csv` and the QC summaries. A sample manifest is generated during fresh-start preprocessing.

## Rheumatoid arthritis: synovial tissue

| Modality | GEO accession | Processed source used |
|---|---|---|
| RNA-seq | GSE112656 | Donor-level count matrix `GSE112656_Counts_matix.txt.gz` |
| ATAC-seq | GSE112654 | Processed donor-level peak files linked from GEO SOFT metadata |
| ChIP-seq | GSE112655 | Processed H3K27ac and H3K4me3 peak files linked from GEO SOFT metadata |

The selected RA cohort contains 10 RA and 10 OA donors with complete profiles for RNA, ATAC, H3K27ac and H3K4me3. The peak-to-gene promoter window uses hg19 annotation with 2,000 bp upstream and 500 bp downstream of the annotated TSS. WGBS, FASTQ/SRA and the large full GEO TAR archive are excluded.

## Acquisition and provenance controls

- `src/scripts/python/00_download_tcga_xena.py` acquires TCGA/Xena data and writes per-cohort download manifests.
- `src/scripts/python/02_download_ra_geo.py` obtains the GEO RNA matrix and parses GEO metadata to retrieve selected processed ATAC/ChIP peak files.
- Download manifests record source URL, size, checksum and acquisition status. `results/QC/provenance_aggregated.csv` contains the aggregated acquisition records from the archived run; local absolute paths have been removed from the submitted copy.
- Raw input files and large processed matrices are not included in this repository. This follows the project's compact-repository policy; the scripts reconstruct them from public sources.
- The public challenge asks entrants to use open or legitimately accessible data and to document sources, selection criteria, labels and processed-file provenance. See the official challenge page: <https://2026.genedis.org/hackathon/>.

## Licensing and access

The pipeline accesses public TCGA/Xena and GEO resources. Reusers are responsible for following the terms and citation requirements of the original repositories and source studies. This repository does not redistribute the large source datasets.
