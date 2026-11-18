#!/usr/bin/env python3
"""Prepare LUAD/LUSC Xena matrices while streaming large HM450 methylation files.

This version preserves the feature-by-sample TSV.GZ output expected by the
MOSAIC-3 downstream R scripts, but does not load the full methylation matrix
into RAM. RNA and CNV continue to use the original in-memory workflow because
they are substantially smaller.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

TCGA_SAMPLE_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z]?)", re.I)
TCGA_CASE_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", re.I)
INVALID_FEATURE_IDS = {"", "nan", "na", "id_ref", "_gene_"}


def locate_layer(raw: Path, layer: str, cohort: str) -> Path:
    candidates = sorted((raw / layer).glob("*.tsv*"))
    if not candidates:
        raise FileNotFoundError(f"Downloaded file not found for {cohort}/{layer}")
    return candidates[0]


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, **kwargs)


def norm_sample(x: str) -> str:
    match = TCGA_SAMPLE_RE.search(str(x).upper())
    return match.group(1) if match else ""


def case_id(x: str) -> str:
    match = TCGA_CASE_RE.search(str(x).upper())
    return match.group(1) if match else ""


def sample_type(x: str) -> str:
    match = re.match(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-(\d{2})[A-Z]?", str(x).upper())
    return match.group(2) if match else ""


def select_primary_samples(ids: list[str]) -> dict[str, str]:
    by_case: dict[str, list[str]] = {}
    for raw_id in ids:
        sample = norm_sample(raw_id)
        if not sample or sample_type(sample) != "01":
            continue
        by_case.setdefault(sample[:12], []).append(sample)

    selected: dict[str, str] = {}
    for case, values in by_case.items():
        values = sorted(
            set(values),
            key=lambda value: (0 if re.search(r"-01A(?:-|$)", value) else 1, value),
        )
        selected[case] = values[0]
    return selected


def load_xena_matrix(path: Path, layer: str) -> pd.DataFrame:
    """Load the relatively small RNA/CNV/SNV matrices using the original logic."""
    print(f"Leyendo {layer}: {path.name}", flush=True)
    header = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    if len(header.columns) < 3:
        raise ValueError(f"Empty or non-tabular Xena matrix: {path}")

    dtype = {str(column): np.float32 for column in header.columns[1:]}
    dtype[str(header.columns[0])] = str
    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=True,
        dtype=dtype,
    )
    if df.empty or df.shape[1] < 3:
        raise ValueError(f"Empty or non-tabular Xena matrix: {path}")

    id_column = df.columns[0]
    first_values = df.iloc[: min(30, len(df)), 0].astype(str).tolist()
    sample_rows = sum(bool(norm_sample(value) or case_id(value)) for value in first_values)

    # Standard Xena matrices have features in rows and sample barcodes in columns.
    if sample_rows >= max(2, len(first_values) // 3):
        df = df.set_index(id_column).T
        df.index = df.index.astype(str)
    else:
        df = df.set_index(id_column)

    sample_columns = [column for column in df.columns if norm_sample(column)]
    if len(sample_columns) < 3:
        sample_row_ids = [value for value in df.index.astype(str) if norm_sample(value)]
        if len(sample_row_ids) >= 3:
            df = df.T
            sample_columns = [column for column in df.columns if norm_sample(column)]
    if len(sample_columns) < 3:
        raise ValueError(
            f"TCGA columns were not detected in matrix {layer}. "
            f"Example columns: {list(df.columns[:8])}"
        )

    df = df.loc[:, sample_columns]
    df.index = df.index.astype(str)
    df = df.loc[~df.index.str.lower().isin(INVALID_FEATURE_IDS)]
    df = df.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.groupby(level=0, sort=False).mean()
    df = df.loc[df.notna().any(axis=1)]
    return df


def map_matrix_to_rna_cases(
    matrix: pd.DataFrame,
    rna_case_map: dict[str, str],
    layer: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    column_map: dict[str, list[str]] = {}
    for column in matrix.columns:
        sample = norm_sample(column)
        case = sample[:12] if sample else ""
        if case in rna_case_map and sample_type(sample) == "01":
            column_map.setdefault(case, []).append(sample)

    chosen: dict[str, str] = {}
    for case, values in column_map.items():
        values = sorted(
            set(values),
            key=lambda value: (0 if re.search(r"-01A(?:-|$)", value) else 1, value),
        )
        chosen[case] = values[0]

    common_cases = [case for case in rna_case_map if case in chosen]
    source_columns: list[str] = []
    for case in common_cases:
        source = chosen[case]
        source_columns.append(next(c for c in matrix.columns if norm_sample(c) == source))

    output = matrix.loc[:, source_columns].copy()
    output.columns = [rna_case_map[case] for case in common_cases]
    report = pd.DataFrame(
        {
            "case_id": list(rna_case_map),
            "sample_id": [rna_case_map[case] for case in rna_case_map],
            f"{layer}_present": [case in chosen for case in rna_case_map],
            f"{layer}_source_sample": [chosen.get(case, "") for case in rna_case_map],
        }
    )
    return output, report


def find_id_col(df: pd.DataFrame) -> str:
    for column in df.columns:
        if column.lower() in {
            "sample", "sample_id", "submitter_id", "bcr_patient_barcode", "_patient"
        }:
            return column
    return df.columns[0]


def first_nonmissing(row: pd.Series, names: list[str]):
    lower = {str(column).lower(): column for column in row.index}
    for name in names:
        column = lower.get(name.lower())
        if column is not None:
            value = row[column]
            if pd.notna(value) and str(value).strip().lower() not in {
                "", "na", "nan", "[not available]", "[not applicable]"
            }:
                return value
    return np.nan


def prepare_clinical(
    raw: Path, rna_case_map: dict[str, str], processed: Path, cohort: str
) -> pd.DataFrame:
    clinical_path = locate_layer(raw, "CLINICAL", cohort)
    clinical = read_table(clinical_path)
    id_column = find_id_col(clinical)
    clinical["__id"] = clinical[id_column].astype(str)
    clinical["case_id"] = clinical["__id"].map(case_id)

    survival_files = list((raw / "CLINICAL").glob(f"*{cohort}_survival*.tsv*"))
    survival = None
    if survival_files:
        survival_table = read_table(survival_files[0])
        survival_id = find_id_col(survival_table)
        survival_table["case_id"] = survival_table[survival_id].astype(str).map(case_id)
        survival = survival_table.drop_duplicates("case_id", keep="first").set_index("case_id")

    clinical = clinical[clinical["case_id"].isin(rna_case_map)].copy()
    clinical = clinical.drop_duplicates("case_id", keep="first").set_index("case_id")
    rows: list[dict[str, Any]] = []

    for current_case, sample in rna_case_map.items():
        row = clinical.loc[current_case] if current_case in clinical.index else pd.Series(dtype=object)
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        survival_row = (
            survival.loc[current_case]
            if survival is not None and current_case in survival.index
            else pd.Series(dtype=object)
        )
        if isinstance(survival_row, pd.DataFrame):
            survival_row = survival_row.iloc[0]

        os_time = first_nonmissing(survival_row, ["OS.time"]) if len(survival_row) else np.nan
        os_event = first_nonmissing(survival_row, ["OS"]) if len(survival_row) else np.nan
        if pd.isna(os_time):
            os_time = first_nonmissing(row, ["OS.time"])
        if pd.isna(os_event):
            os_event = first_nonmissing(row, ["OS"])

        if isinstance(os_event, str):
            event_text = os_event.strip().lower()
            event = (
                1 if event_text in {"dead", "deceased", "1", "yes", "true"} or "dead" in event_text
                else 0 if event_text in {"alive", "0", "no", "false"} or "alive" in event_text
                else np.nan
            )
        else:
            try:
                event_value = float(os_event) if pd.notna(os_event) else np.nan
                event = int(event_value) if np.isfinite(event_value) and event_value in (0, 1) else np.nan
            except (TypeError, ValueError):
                event = np.nan

        try:
            time_days = float(os_time)
        except (TypeError, ValueError):
            time_days = np.nan

        # Fallback to follow-up fields; days_to_death is not itself an event label.
        days_to_death = first_nonmissing(row, ["days_to_death"])
        days_last_followup = first_nonmissing(row, ["days_to_last_followup", "days_to_last_follow_up"])
        vital_status = first_nonmissing(row, ["vital_status"])
        if pd.isna(event) and pd.notna(vital_status):
            vital_text = str(vital_status).strip().lower()
            event = 1 if vital_text in {"dead", "deceased"} else 0 if vital_text == "alive" else np.nan

        try:
            if (not np.isfinite(time_days) or time_days <= 0) and pd.notna(days_to_death) and float(days_to_death) > 0:
                time_days = float(days_to_death)
                event = 1 if pd.isna(event) else event
            elif (not np.isfinite(time_days) or time_days <= 0) and pd.notna(days_last_followup) and float(days_last_followup) > 0:
                time_days = float(days_last_followup)
                event = 0 if pd.isna(event) else event
        except (TypeError, ValueError):
            pass

        age = first_nonmissing(row, ["age_at_initial_pathologic_diagnosis", "age_at_diagnosis", "age"])
        try:
            age = float(age)
        except (TypeError, ValueError):
            age = np.nan
        if pd.notna(age) and age > 200:
            age /= 365.25

        sex = first_nonmissing(row, ["gender", "sex"])
        stage = first_nonmissing(row, ["pathologic_stage", "ajcc_pathologic_tumor_stage", "tumor_stage"])
        rows.append(
            {
                "case_id": current_case,
                "sample_id": sample,
                "time_days": time_days,
                "event": event,
                "age_years_candidate": age,
                "sex_candidate": sex,
                "stage_group": stage,
            }
        )

    output = pd.DataFrame(rows).replace([np.inf, -np.inf], np.nan)
    output.to_csv(processed / "clinical_all.csv", index=False)
    return output


def prepare_snv(raw: Path, canonical_cases: dict[str, str]) -> pd.DataFrame | None:
    files = list((raw / "SNV").glob("*.tsv*"))
    if not files:
        return None
    raw_matrix = load_xena_matrix(files[0], "SNV_gene_level")
    matrix, _ = map_matrix_to_rna_cases(raw_matrix, canonical_cases, "SNV")
    return matrix


def _choose_methylation_columns(
    header: list[str],
    canonical_cases: dict[str, str],
    eligible_cases: list[str],
) -> tuple[list[int], list[str], list[str], pd.DataFrame]:
    """Choose one primary-tumour methylation column per eligible RNA case."""
    wanted = set(eligible_cases)
    candidates: dict[str, list[tuple[str, int, str]]] = {}
    for column_index, column_name in enumerate(header[1:], start=1):
        sample = norm_sample(column_name)
        if not sample or sample_type(sample) != "01":
            continue
        current_case = sample[:12]
        if current_case in wanted:
            candidates.setdefault(current_case, []).append((sample, column_index, column_name))

    selected_cases: list[str] = []
    selected_indices: list[int] = []
    selected_output_samples: list[str] = []
    report_rows: list[dict[str, Any]] = []

    for current_case in eligible_cases:
        options = candidates.get(current_case, [])
        if not options:
            report_rows.append(
                {
                    "case_id": current_case,
                    "sample_id": canonical_cases[current_case],
                    "METHYLATION_present": False,
                    "METHYLATION_source_sample": "",
                }
            )
            continue
        options.sort(key=lambda item: (0 if re.search(r"-01A(?:-|$)", item[0]) else 1, item[0]))
        chosen_sample, chosen_index, _ = options[0]
        selected_cases.append(current_case)
        selected_indices.append(chosen_index)
        selected_output_samples.append(canonical_cases[current_case])
        report_rows.append(
            {
                "case_id": current_case,
                "sample_id": canonical_cases[current_case],
                "METHYLATION_present": True,
                "METHYLATION_source_sample": chosen_sample,
            }
        )

    return (
        selected_cases,
        selected_indices,
        selected_output_samples,
        pd.DataFrame(report_rows),
    )


def stream_methylation_matrix(
    path: Path,
    output_path: Path,
    canonical_cases: dict[str, str],
    eligible_cases: list[str],
    chunk_rows: int = 5000,
) -> tuple[list[str], int, pd.DataFrame]:
    """Stream HM450 rows to a compressed TSV without loading the matrix into RAM.

    Pandas only holds one small chunk (default 5,000 probes) at a time. This
    preserves the feature-by-sample layout expected by the downstream R code.
    """
    print(f"Streaming METHYLATION in chunks of {chunk_rows:,} rows: {path.name}", flush=True)
    header_frame = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    header = [str(column) for column in header_frame.columns]
    if len(header) < 3:
        raise ValueError(f"Non-tabular methylation matrix: {path}")

    selected_cases, selected_indices, output_samples, report = _choose_methylation_columns(
        header, canonical_cases, eligible_cases
    )
    if len(selected_cases) < 3:
        raise ValueError(
            f"Only {len(selected_cases)} eligible primary-tumour methylation samples were found in {path}"
        )
    if len(set(selected_indices)) != len(selected_indices):
        raise ValueError("Methylation column selection produced duplicated source columns")

    id_column = header[0]
    selected_source_columns = [header[index] for index in selected_indices]
    use_columns = [id_column, *selected_source_columns]
    dtype_map = {id_column: str}
    dtype_map.update({column: np.float32 for column in selected_source_columns})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(output_path.name + ".partial")
    seen_features: set[str] = set()
    feature_count = 0
    skipped_empty = 0
    duplicate_examples: list[str] = []
    wrote_header = False

    try:
        with gzip.open(temporary_path, "wt", encoding="utf-8", newline="", compresslevel=1) as out_handle:
            reader = pd.read_csv(
                path,
                sep="\t",
                compression="infer",
                usecols=use_columns,
                dtype=dtype_map,
                chunksize=chunk_rows,
                low_memory=True,
            )
            for chunk_number, chunk in enumerate(reader, start=1):
                # read_csv(usecols=...) may not preserve the order supplied in usecols.
                chunk = chunk.loc[:, use_columns]
                feature_ids = chunk[id_column].fillna("").astype(str).str.strip()
                values = chunk.loc[:, selected_source_columns].to_numpy(dtype=np.float32, copy=False)
                values[~np.isfinite(values)] = np.nan
                row_has_numeric = np.isfinite(values).any(axis=1)

                keep = np.zeros(len(chunk), dtype=bool)
                for row_index, feature_id in enumerate(feature_ids.tolist()):
                    if feature_id.lower() in INVALID_FEATURE_IDS or not row_has_numeric[row_index]:
                        if not row_has_numeric[row_index]:
                            skipped_empty += 1
                        continue
                    if feature_id in seen_features:
                        if len(duplicate_examples) < 10:
                            duplicate_examples.append(feature_id)
                        continue
                    seen_features.add(feature_id)
                    keep[row_index] = True

                if not keep.any():
                    continue

                output_chunk = chunk.loc[keep, use_columns].copy()
                output_chunk.columns = ["feature_id", *output_samples]
                output_chunk.to_csv(
                    out_handle,
                    sep="\t",
                    index=False,
                    header=not wrote_header,
                    float_format="%.9g",
                    na_rep="NA",
                    lineterminator="\n",
                )
                wrote_header = True
                feature_count += int(keep.sum())
                if chunk_number % 10 == 0:
                    print(
                        f"  processed approximately {chunk_number * chunk_rows:,} probes; "
                        f"written {feature_count:,}",
                        flush=True,
                    )

        if duplicate_examples:
            raise ValueError(
                "Duplicate methylation feature IDs were detected (examples: "
                + ", ".join(duplicate_examples)
                + "). The streaming writer stopped rather than silently changing the original mean-aggregation behavior."
            )
        if feature_count == 0 or not wrote_header:
            raise ValueError(f"No usable methylation features were written from {path}")
        os.replace(temporary_path, output_path)
    except Exception:
        try:
            temporary_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise

    print(
        f"METHYLATION streamed: {feature_count:,} features; "
        f"{len(output_samples)} aligned samples; {skipped_empty:,} all-missing rows skipped.",
        flush=True,
    )
    return output_samples, feature_count, report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--cohort", choices=["LUAD", "LUSC"], required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    cohort = args.cohort
    raw = root / "data" / "raw" / cohort
    processed = root / "data" / "processed" / cohort
    processed.mkdir(parents=True, exist_ok=True)

    # RNA defines the canonical primary-tumour sample per case.
    rna_raw = load_xena_matrix(locate_layer(raw, "RNA", cohort), "RNA")
    rna_case_map = select_primary_samples(list(rna_raw.columns))
    if len(rna_case_map) < 100:
        raise RuntimeError(
            f"{cohort}: only {len(rna_case_map)} primary cases detected in RNA; check identifiers."
        )
    print(f"{cohort} RNA primary tumour: {len(rna_case_map)} casos.", flush=True)

    rna_source_columns = [
        next(column for column in rna_raw.columns if norm_sample(column) == sample)
        for sample in rna_case_map.values()
    ]
    rna = rna_raw.loc[:, rna_source_columns].copy()
    rna.columns = list(rna_case_map.values())
    del rna_raw

    # CNV is smaller and can be loaded normally.
    cnv_raw = load_xena_matrix(locate_layer(raw, "CNV", cohort), "CNV")
    cnv, cnv_report = map_matrix_to_rna_cases(cnv_raw, rna_case_map, "CNV")
    del cnv_raw
    print(f"CNV: {cnv.shape[1]} pacientes alineados, {cnv.shape[0]} características.", flush=True)

    # Resolve valid survival cases before streaming methylation, so only relevant
    # primary-tumour columns are read from the very large HM450 file.
    clinical = prepare_clinical(raw, rna_case_map, processed, cohort)
    valid_samples = set(
        clinical.loc[
            clinical.time_days.notna()
            & clinical.event.notna()
            & (clinical.time_days > 0),
            "sample_id",
        ].astype(str)
    )
    rna_samples = set(rna.columns.astype(str))
    cnv_samples = set(cnv.columns.astype(str))
    eligible_cases = [
        case
        for case, sample in rna_case_map.items()
        if sample in valid_samples and sample in rna_samples and sample in cnv_samples
    ]

    methylation_path = locate_layer(raw, "METHYLATION", cohort)
    methylation_output = processed / "METHYLATION.tsv.gz"
    methylation_samples, methylation_features, methylation_report = stream_methylation_matrix(
        methylation_path,
        methylation_output,
        rna_case_map,
        eligible_cases,
    )

    # Streaming output columns are already ordered and renamed to the canonical
    # RNA sample IDs. Their order follows eligible_cases / the RNA cohort order.
    common = methylation_samples
    if len(common) < 100:
        raise RuntimeError(
            f"{cohort}: only {len(common)} complete cases across RNA, CNV, methylation and survival."
        )

    rna = rna.loc[:, common]
    cnv = cnv.loc[:, common]
    clinical = clinical[clinical.sample_id.isin(common)].copy().sort_values("sample_id")

    # Optional somatic SNV block. It is not part of the primary survival benchmark.
    snv_error = ""
    snv_features = 0
    try:
        snv = prepare_snv(raw, rna_case_map)
    except Exception as exc:
        snv = None
        snv_error = str(exc)
        print(f"SNV exploratorio omitido (no bloquea el benchmark principal): {exc}", flush=True)
    if snv is not None:
        snv.to_csv(
            processed / "SNV_exploratory_candidates.tsv.gz",
            sep="\t",
            compression={"method": "gzip", "compresslevel": 1},
            index_label="feature_id",
        )
        snv_features = snv.shape[0]
        print(
            f"SNV exploratorio: {snv.shape[1]} pacientes; excluido del benchmark primario "
            "hasta verificar callability.",
            flush=True,
        )
        del snv

    pd.DataFrame(
        [
            {
                "status": "prepared" if snv_features > 0 else ("unavailable_or_skipped" if not snv_error else "parse_failed"),
                "features": snv_features,
                "error": snv_error,
                "interpretation": "Exploratory only; missing calls are not interpreted as wild type without callability.",
            }
        ]
    ).to_csv(processed / "SNV_exploratory_QC.csv", index=False)

    # Save compact matrices. HM450 was already streamed to its final output file.
    rna.to_csv(
        processed / "RNA.tsv.gz",
        sep="\t",
        compression={"method": "gzip", "compresslevel": 1},
        index_label="feature_id",
    )
    cnv.to_csv(
        processed / "CNV.tsv.gz",
        sep="\t",
        compression={"method": "gzip", "compresslevel": 1},
        index_label="feature_id",
    )
    clinical.to_csv(processed / "clinical.csv", index=False)

    sample_qc = pd.DataFrame({"sample_id": common, "case_id": [sample[:12] for sample in common]})
    sample_qc["RNA_present"] = True
    sample_qc["CNV_present"] = True
    sample_qc["METHYLATION_present"] = True
    sample_qc.to_csv(processed / "sample_alignment.csv", index=False)
    clinical[["case_id", "sample_id", "time_days", "event"]].sort_values(
        "sample_id"
    ).to_csv(processed / "selected_sample_manifest.csv", index=False)

    pd.concat([cnv_report, methylation_report], ignore_index=True, sort=False).to_csv(
        processed / "layer_alignment_report.csv", index=False
    )

    pd.DataFrame(
        [
            {
                "cohort": f"TCGA-{cohort}",
                "n_primary_RNA_cases": len(rna_case_map),
                "n_model_complete_cases": len(common),
                "n_events": int(clinical.event.sum()),
                "n_censored": int((clinical.event == 0).sum()),
                "RNA_features": rna.shape[0],
                "CNV_features": cnv.shape[0],
                "METHYLATION_features": methylation_features,
                "SNV_exploratory_features": snv_features,
            }
        ]
    ).to_csv(processed / "cohort_qc.csv", index=False)

    print(f"\n{cohort} prepared. Complete cases across all three layers and survival: {len(common)}", flush=True)
    print(f"Eventos: {int(clinical.event.sum())}; censurados: {int((clinical.event == 0).sum())}", flush=True)
    print(
        f"Se conservan RNA={rna.shape}, CNV={cnv.shape}, "
        f"Methylation=({methylation_features}, {len(common)}).",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
