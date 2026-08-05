#!/usr/bin/env python3
"""Normalize the existing PRSS55 window-prediction table for structural analysis."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED = {"seq_id", "protein_name", "window_idx", "start_pos", "end_pos", "8mer_sequence", "prob_cleavage", "prediction"}


def accession(seq_id):
    fields = seq_id.split("|")
    if len(fields) == 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    raise ValueError(f"Cannot extract UniProt accession from seq_id: {seq_id}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path, help="TSV or CSV in the supplied prediction-table format")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--p1-offset", type=int, choices=range(1, 8),
                    help="Known P1 position within each 8-mer; omit only when the cleavage bond is unresolved")
    args = ap.parse_args()
    delimiter = "\t" if args.input.suffix.lower() in {".tsv", ".tab"} else ","
    with args.input.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or not REQUIRED.issubset(reader.fieldnames):
            raise ValueError("Input is missing: " + ", ".join(sorted(REQUIRED - set(reader.fieldnames or []))))
        normalized = []
        for line, row in enumerate(reader, start=2):
            sequence = row["8mer_sequence"].strip().upper()
            start, end = int(row["start_pos"]), int(row["end_pos"])
            if len(sequence) != 8 or end - start + 1 != 8:
                raise ValueError(f"Line {line}: positions and 8-mer length disagree")
            normalized.append({"accession": accession(row["seq_id"].strip()), "site_start": start,
                               "eight_mer": sequence, "p1_offset": args.p1_offset or "", "model_score": row["prob_cleavage"],
                               "source_isoform": "", "source_id": f"{row['seq_id']}|window_{row['window_idx']}"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "site_start", "eight_mer", "p1_offset", "model_score", "source_isoform", "source_id"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(normalized)
    print(f"Wrote {len(normalized)} normalized predictions to {args.output}")


if __name__ == "__main__":
    main()
