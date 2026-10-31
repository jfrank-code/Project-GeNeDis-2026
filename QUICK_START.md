# Quick start

Requirements: Python 3.10+, R 4.3+, internet connection, and sufficient disk space/RAM.

From PowerShell in this directory:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\run_mosaic3.ps1
```

The pipeline automatically downloads LUAD/LUSC from UCSC Xena and RA data, including processed ATAC and ChIP peak files, from GEO. No manual ATAC folder copy is required. See `README.md` and `docs/REPRODUCIBILITY.md` for the complete workflow.


Optional PDF rendering: install `requirements_pdf_optional.txt` to regenerate the PDF using WeasyPrint. If unavailable, the report generator falls back to ReportLab; the included PDF was rendered from the HTML report with WeasyPrint.
