# GitHub upload checklist

## Before uploading

- Upload the contents of this project folder (or extract the submission ZIP first).
- Keep `data/raw/` and `data/processed/` empty except for their `.gitkeep` files. The `.gitignore` already excludes downloaded inputs and large regenerated matrices.
- Keep `results/` and `Reports/`: these contain the archived benchmark outputs, QC evidence, and final English report.
- Do not add API keys, credentials, raw subject-level source files, or local absolute-path logs.

## Upload through the GitHub website

1. Create a new public repository named `MOSAIC-3` (or another chosen name).
2. Extract `MOSAIC3_GeNeDis_2026_Submission.zip` locally.
3. Upload the project files and folders from the extracted directory to the repository root.
4. Confirm that `README.md`, `run_mosaic3.ps1`, `src/`, `results/`, `Reports/`, `docs/`, and `config/` appear at the repository root.
5. Open `Reports/MOSAIC3_final_report.pdf` from GitHub and verify that it renders correctly.

## Upload with Git (optional)

```bash
git init
git add .
git commit -m "Submit MOSAIC-3 GeNeDis 2026"
git branch -M main
git remote add origin https://github.com/<YOUR-ACCOUNT>/<YOUR-REPOSITORY>.git
git push -u origin main
```

Replace the placeholder remote URL with the repository URL created in GitHub. No data download is required before publishing the code and archived results.
