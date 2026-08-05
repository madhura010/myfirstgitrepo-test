#!/usr/bin/env python3
"""Extract auditable, site-level features for predicted PRSS55 cleavage sites.

This script deliberately does not label sites as substrates or derive a cleavage
probability. It validates the input against UniProt, retrieves annotations, and
optionally computes AlphaFold-model-derived residue RSA using the installed
aganitha-uniprot-parser implementation.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
PARSER_PATH = ROOT / "skills" / "aganitha-uniprot-parser" / "scripts" / "fetch_uniprot.py"
OUTPUT_COLUMNS = (
    "source_id accession resolved_accession source_isoform site_start site_end eight_mer "
    "p1_position p1prime_position model_score sequence_verified protein_name gene_name "
    "processing_overlap mature_peptide_overlap chain_annotations transmembrane_overlap topology_annotations domain_annotations family_annotations "
    "glycosylation_count_near_site modified_residue_count_near_site disulfide_overlap rsa_p1 rsa_p1prime rsa_window_mean "
    "rsa_window_min plddt_window_mean ss_p1 ss_p1prime ss_window_pattern ss_helix_fraction ss_strand_fraction ss_loop_fraction "
    "structure_source evidence_gaps warnings retrieved_at"
).split()


def load_parser():
    if not PARSER_PATH.exists():
        raise RuntimeError(f"Required parser not found: {PARSER_PATH}")
    spec = importlib.util.spec_from_file_location("uniprot_parser", PARSER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def read_sites(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"accession", "site_start", "eight_mer", "model_score"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Input requires columns: " + ", ".join(sorted(required)))
        for line, row in enumerate(reader, start=2):
            if not row.get("accession", "").strip() or row["accession"].lstrip().startswith("#"):
                continue
            try:
                start = int(row["site_start"])
                score = float(row["model_score"])
            except ValueError as exc:
                raise ValueError(f"Line {line}: site_start and model_score must be numeric") from exc
            offset_text = row.get("p1_offset", "").strip()
            try:
                offset = int(offset_text) if offset_text else None
            except ValueError as exc:
                raise ValueError(f"Line {line}: p1_offset must be blank or an integer 1..7") from exc
            peptide = row["eight_mer"].strip().upper()
            if start < 1 or (offset is not None and offset not in range(1, 8)) or len(peptide) != 8 or not peptide.isalpha():
                raise ValueError(f"Line {line}: require 1-based start, optional p1_offset 1..7, and an 8-aa sequence")
            yield line, {**row, "site_start": start, "p1_offset": offset, "model_score": score, "eight_mer": peptide}


def bounds(feature, parser):
    return parser.feature_bounds(feature)


def overlaps(start, end, feature, parser):
    left, right = bounds(feature, parser)
    return left is not None and left <= end and right >= start


def descriptions(features):
    """Keep the feature type when UniProt intentionally omits a description.

    A blank description for a Signal or Propeptide feature is still meaningful
    and must not be serialised as the same sentinel used for 'no overlap'.
    """
    labels = set()
    for feature in features:
        kind = feature.get("type", "Feature")
        detail = feature.get("description", "") or ""
        labels.add(f"{kind}: {detail}" if detail else kind)
    return "; ".join(sorted(labels)) or "-"


def plddt_by_residue(pdb_text: str):
    """Read AlphaFold pLDDT stored in CA-atom B factors; returns residue -> pLDDT."""
    values = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            values[int(line[22:26])] = float(line[60:66])
        except ValueError:
            continue
    return values


def secondary_structure_by_residue(pdb_text: str):
    """Assign C3 DSSP states from backbone coordinates with local PyDSSP.

    The output uses H (alpha helix), E (beta strand), and C (loop/other). This
    is a structure assignment for the selected AlphaFold model, not an
    experimentally observed ensemble or a sequence-only prediction.
    """
    vendor = ROOT / ".vendor"
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    try:
        import pydssp
    except ImportError as exc:
        raise RuntimeError("PyDSSP unavailable; install project dependencies") from exc
    try:
        coordinates = pydssp.read_pdbtext(pdb_text)
        states = pydssp.assign(coordinates, out_type="c3")
    except (AssertionError, ValueError, IndexError) as exc:
        raise RuntimeError(f"PyDSSP assignment failed: {exc}") from exc
    residue_numbers, seen = [], set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "N":
            continue
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        if number not in seen:
            residue_numbers.append(number); seen.add(number)
    if len(residue_numbers) != len(states):
        raise RuntimeError("PyDSSP residue count does not match PDB residue numbering")
    conversion = {"H": "H", "E": "E", "-": "C"}
    return {number: conversion.get(str(state), "C") for number, state in zip(residue_numbers, states)}


def family_annotations(entry):
    """Return curated InterPro/Pfam memberships without inventing coordinates."""
    labels = set()
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") not in {"InterPro", "Pfam"}:
            continue
        name = next((p.get("value") for p in ref.get("properties", []) if p.get("key") == "EntryName"), "")
        labels.add(f"{ref['database']}: {name or ref.get('id', '-')}")
    return "; ".join(sorted(labels)) or "-"


def cached_entry(accession, parser, cache_dir: Path):
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{accession.upper()}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    entry = parser.entry_for(accession)
    cache_path.write_text(json.dumps(entry, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return entry


def structural_features(accession, parser):
    """Return RSA/pLDDT evidence. Failures are reported, never converted to zero exposure."""
    try:
        rsa, pdb_url = parser.alphafold_rsa(accession)
        pdb_text = parser.get(pdb_url, "text/plain")
        plddt = plddt_by_residue(pdb_text)
        secondary = secondary_structure_by_residue(pdb_text)
        return rsa, plddt, secondary, pdb_url, None
    except RuntimeError as exc:
        return {}, {}, {}, "", str(exc)


def feature_row(row, parser, cache_dir, with_rsa):
    entry = cached_entry(row["accession"], parser, cache_dir)
    accession = entry["primaryAccession"]
    sequence = entry.get("sequence", {}).get("value", "")
    start, end = row["site_start"], row["site_start"] + 7
    offset = row["p1_offset"]
    p1 = start + offset - 1 if offset is not None else None
    p1prime = start + offset if offset is not None else None
    warnings, gaps = [], []
    observed = sequence[start - 1:end] if end <= len(sequence) else ""
    verified = observed == row["eight_mer"]
    if not verified:
        warnings.append(f"sequence mismatch (UniProt has '{observed or 'out of range'}')")

    features = entry.get("features", [])
    by_type = lambda kinds: [f for f in features if f.get("type") in kinds and overlaps(start, end, f, parser)]
    # A UniProt `Chain` commonly spans the mature protein and must not be treated
    # as evidence that a site is removed by processing. Only features that can
    # remove or release the local sequence are used as a processing flag.
    processing = by_type({"Signal", "Propeptide", "Initiator methionine"})
    # `Peptide` is an annotated mature cleavage product. Its overlap supports
    # the existence of the local sequence; it does not mean the sequence was
    # removed and must never be used as an exclusion criterion.
    mature_peptide = by_type({"Peptide"})
    chain = by_type({"Chain"})
    membrane = by_type({"Transmembrane", "Intramembrane"})
    topology = by_type({"Topological domain"})
    domain = by_type({"Domain", "Region", "Repeat", "Motif", "Zinc finger", "Coiled coil"})
    disulfide = by_type({"Disulfide bond"})
    radius = 8
    glyco, modified = [], []
    for f in features:
        if f.get("type") not in {"Glycosylation", "Modified residue"}:
            continue
        left, right = bounds(f, parser)
        if left is not None and left <= end + radius and right >= start - radius:
            (glyco if f.get("type") == "Glycosylation" else modified).append(f)

    rsa, plddt, secondary, structure_source, structure_error = ({}, {}, {}, "", None)
    if with_rsa:
        rsa, plddt, secondary, structure_source, structure_error = structural_features(accession, parser)
        if structure_error:
            gaps.append("RSA/pLDDT unavailable: " + structure_error)
    else:
        gaps.append("RSA/pLDDT not requested; rerun with --with-rsa")
    window_rsa = [rsa[p] for p in range(start, end + 1) if p in rsa]
    window_plddt = [plddt[p] for p in range(start, end + 1) if p in plddt]
    window_ss = [secondary[p] for p in range(start, end + 1) if p in secondary]
    if with_rsa and len(window_rsa) != 8:
        gaps.append("incomplete structural coverage for 8-mer")
    if with_rsa and len(window_ss) != 8:
        gaps.append("incomplete secondary-structure coverage for 8-mer")
    protein = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "-")
    genes = entry.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "-") if genes else "-"
    return {
        "source_id": row.get("source_id", "") or f"line_{row.get('_line', '')}",
        "accession": row["accession"], "resolved_accession": accession,
        "source_isoform": row.get("source_isoform", ""), "site_start": start, "site_end": end,
        "eight_mer": row["eight_mer"], "p1_position": p1 or "", "p1prime_position": p1prime or "",
        "model_score": row["model_score"], "sequence_verified": str(verified).lower(),
        "protein_name": protein, "gene_name": gene, "processing_overlap": descriptions(processing),
        "mature_peptide_overlap": descriptions(mature_peptide),
        "chain_annotations": descriptions(chain),
        "transmembrane_overlap": descriptions(membrane), "topology_annotations": descriptions(topology),
        "domain_annotations": descriptions(domain), "family_annotations": family_annotations(entry),
        "glycosylation_count_near_site": len(glyco),
        "modified_residue_count_near_site": len(modified),
        "disulfide_overlap": str(bool(disulfide)).lower(), "rsa_p1": rsa.get(p1, "") if p1 else "",
        "rsa_p1prime": rsa.get(p1prime, "") if p1prime else "",
        "rsa_window_mean": f"{mean(window_rsa):.3f}" if window_rsa else "",
        "rsa_window_min": f"{min(window_rsa):.3f}" if window_rsa else "",
        "plddt_window_mean": f"{mean(window_plddt):.3f}" if window_plddt else "",
        "ss_p1": secondary.get(p1, "") if p1 else "", "ss_p1prime": secondary.get(p1prime, "") if p1prime else "",
        "ss_window_pattern": "".join(window_ss) if len(window_ss) == 8 else "",
        "ss_helix_fraction": f"{window_ss.count('H') / len(window_ss):.3f}" if window_ss else "",
        "ss_strand_fraction": f"{window_ss.count('E') / len(window_ss):.3f}" if window_ss else "",
        "ss_loop_fraction": f"{window_ss.count('C') / len(window_ss):.3f}" if window_ss else "",
        "structure_source": structure_source, "evidence_gaps": " | ".join(gaps) or "-",
        "warnings": " | ".join(warnings) or "-", "retrieved_at": dt.date.today().isoformat(),
    }


def main():
    parser_args = argparse.ArgumentParser(description=__doc__)
    parser_args.add_argument("--input", required=True, type=Path, help="Predicted-site CSV; use data/prss55_predicted_sites.template.csv")
    parser_args.add_argument("--output", required=True, type=Path, help="Output TSV")
    parser_args.add_argument("--cache-dir", type=Path, default=ROOT / "cache" / "uniprot")
    parser_args.add_argument("--with-rsa", action="store_true", help="Fetch AlphaFold models and compute RSA/pLDDT; slower")
    args = parser_args.parse_args()
    parser = load_parser()
    rows = []
    for line, site in read_sites(args.input):
        site["_line"] = line
        try:
            rows.append(feature_row(site, parser, args.cache_dir, args.with_rsa))
        except RuntimeError as exc:
            rows.append({column: "" for column in OUTPUT_COLUMNS} | {"source_id": site.get("source_id", f"line_{line}"), "accession": site["accession"], "warnings": str(exc), "evidence_gaps": "UniProt retrieval failed"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} site rows to {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
