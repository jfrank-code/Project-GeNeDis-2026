[CmdletBinding()]
param(
    [switch]$ForceDownloads,
    [switch]$SkipRDependencyInstall
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Root

function Invoke-Step {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Label,

        [Parameter(Mandatory = $true)]
        [string]$Exe,

        [Parameter(Mandatory = $true)]
        [string[]]$StepArguments
    )

    Write-Host "`n============================================================" -ForegroundColor DarkCyan
    Write-Host $Label -ForegroundColor Cyan
    Write-Host ("{0} {1}" -f $Exe, ($StepArguments -join ' ')) -ForegroundColor DarkGray

    & $Exe @StepArguments
    $exitCode = $LASTEXITCODE

    if ($exitCode -ne 0) {
        throw "Step '$Label' failed with exit code $exitCode. Completed downloads are preserved; fix the issue and rerun the same command."
    }
}

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    throw "Python 3.10+ was not found on PATH."
}
if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
    throw "Rscript was not found. Install R 4.3+ and add Rscript to PATH."
}

Write-Host "MOSAIC-3 | FRESH START | Root: $Root" -ForegroundColor Green
Write-Host "No pre-existing local matrices or prior pipeline installation are required." -ForegroundColor Green

Invoke-Step -Label "Install Python dependencies" -Exe "python" -StepArguments @("-m", "pip", "install", "-r", "requirements.txt")

if (-not $SkipRDependencyInstall) {
    Invoke-Step -Label "Install/check R and Bioconductor dependencies" -Exe "Rscript" -StepArguments @("src/scripts/R/00_install_dependencies.R")
}

$tcgaArgs = @("src/scripts/python/00_download_tcga_xena.py", "--root", ".")
if ($ForceDownloads) {
    $tcgaArgs += "--force"
}
Invoke-Step -Label "Download LUAD + LUSC from UCSC Xena: RNA, CNV, HM450, clinical and survival data" -Exe "python" -StepArguments $tcgaArgs

foreach ($cohort in @("LUAD", "LUSC")) {
    Invoke-Step -Label "Prepare and align TCGA-${cohort}: primary tumour, one sample per case, three layers + OS" -Exe "python" -StepArguments @("src/scripts/python/03_prepare_tcga_matrices.py", "--root", ".", "--cohort", $cohort)
}

$raArgs = @("src/scripts/python/02_download_ra_geo.py", "--root", ".")
if ($ForceDownloads) {
    $raArgs += "--force"
}
Invoke-Step -Label "Download RA from GEO: RNA-seq, ATAC-seq, H3K27ac and H3K4me3 (excluding WGBS/FASTQ/full TAR)" -Exe "python" -StepArguments $raArgs

Invoke-Step -Label "Prepare donor-level RA profiles and audit the cross-modality intersection" -Exe "Rscript" -StepArguments @("src/scripts/R/01_prepare_ra_features.R")
Invoke-Step -Label "Preflight: provenance, identifiers, omics layers and endpoints" -Exe "python" -StepArguments @("src/scripts/python/05_validate_project.py", "--root", ".")
Invoke-Step -Label "LUAD + LUSC survival benchmark (nested validation)" -Exe "Rscript" -StepArguments @("src/scripts/R/02_run_nsclc_survival_benchmark.R")
Invoke-Step -Label "RA/OA multi-omics benchmark and biomarker stability" -Exe "Rscript" -StepArguments @("src/scripts/R/03_run_ra_benchmark.R")
Invoke-Step -Label "RA-LUAD/LUSC cross-indication comparison and GO enrichment" -Exe "Rscript" -StepArguments @("src/scripts/R/04_cross_indication.R")
Invoke-Step -Label "Generate HTML and PDF reports from archived results" -Exe "python" -StepArguments @("src/scripts/python/04_generate_report.py", "--root", ".")

Write-Host "`nMOSAIC-3 COMPLETED." -ForegroundColor Green
Write-Host "Report:  $Root\Reports\MOSAIC3_final_report.html" -ForegroundColor Yellow
Write-Host "QC:      $Root\results\QC\preflight_summary.json" -ForegroundColor Yellow

