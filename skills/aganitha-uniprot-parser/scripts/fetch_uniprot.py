#!/usr/bin/env python3
"""Retrieve UniProt features and calculate model-derived RSA from AlphaFold PDB files.

Uses only the Python standard library so a fresh Codex environment can run it.
"""
import argparse
import csv
import datetime as dt
import json
import math
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

UNIPROT = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction/"
FEATURES = {"Domain", "Region", "Repeat", "Motif", "Zinc finger", "Coiled coil", "Topological domain", "Transmembrane", "Intramembrane"}
RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
# Maximum solvent-accessible areas (Å²), Tien et al. 2013, complete Gly-X-Gly context.
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
           "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
           "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
           "TYR": 263, "VAL": 174, "SEC": 167, "PYL": 236}
ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")


def get(url, accept="application/json"):
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "aganitha-uniprot-parser/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            return response.read().decode("utf-8")
    except urllib.error.HTTPError as err:
        detail = err.read().decode("utf-8", "replace")[:300]
        raise RuntimeError(f"HTTP {err.code} for {url}: {detail}") from err
    except urllib.error.URLError as err:
        raise RuntimeError(f"Could not reach {url}: {err.reason}") from err


def entry_for(token):
    token = token.strip().upper()
    try:
        return json.loads(get(f"{UNIPROT}/{urllib.parse.quote(token)}.json"))
    except RuntimeError as direct_error:
        query = urllib.parse.urlencode({"query": f"(accession:{token} OR id:{token})", "format": "json", "size": 1})
        data = json.loads(get(f"{UNIPROT}/search?{query}"))
        if not data.get("results"):
            raise RuntimeError(f"Could not resolve '{token}'. {direct_error}")
        return data["results"][0]


def feature_bounds(feature):
    loc = feature.get("location", {})
    try:
        return int(loc["start"]["value"]), int(loc["end"]["value"])
    except (KeyError, TypeError, ValueError):
        return None, None


def selected_features(entry):
    rows = []
    for feature in entry.get("features", []):
        if feature.get("type") not in FEATURES:
            continue
        start, end = feature_bounds(feature)
        if start is None:
            continue
        rows.append({"type": feature["type"], "start": start, "end": end,
                     "description": feature.get("description", "-") or "-",
                     "id": feature.get("featureId", "-") or "-"})
    return rows


def fibonacci_sphere(count):
    points = []
    golden = math.pi * (3 - math.sqrt(5))
    for i in range(count):
        y = 1 - (i / float(count - 1)) * 2
        radius = math.sqrt(max(0, 1 - y * y))
        theta = golden * i
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points


def pdb_atoms(pdb):
    atoms, seen = [], set()
    for line in pdb.splitlines():
        if not line.startswith("ATOM") or (len(line) > 16 and line[16] not in (" ", "A")):
            continue
        try:
            atom = line[12:16].strip()
            residue = line[17:20].strip().upper()
            number = int(line[22:26])
            element = (line[76:78].strip().upper() or atom[0].upper())
            key = (number, residue, atom)
            if key in seen:
                continue
            seen.add(key)
            atoms.append((number, residue, float(line[30:38]), float(line[38:46]), float(line[46:54]), RADII.get(element, 1.70)))
        except (IndexError, ValueError):
            continue
    if not atoms:
        raise RuntimeError("AlphaFold file contained no parsable ATOM records")
    return atoms


def residue_sasa(atoms, points=240, probe=1.4):
    """All-atom Shrake-Rupley SASA with a spatial hash for nearby occluders."""
    cell = 6.6  # > twice the largest expanded atomic radius
    grid = defaultdict(list)
    expanded = [atom[5] + probe for atom in atoms]
    for i, atom in enumerate(atoms):
        grid[(int(atom[2] // cell), int(atom[3] // cell), int(atom[4] // cell))].append(i)
    sphere, output = fibonacci_sphere(points), defaultdict(float)
    for i, atom in enumerate(atoms):
        x, y, z, radius = atom[2], atom[3], atom[4], expanded[i]
        bucket = (int(x // cell), int(y // cell), int(z // cell))
        nearby = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nearby.extend(grid[(bucket[0] + dx, bucket[1] + dy, bucket[2] + dz)])
        available = 0
        for px, py, pz in sphere:
            sx, sy, sz = x + radius * px, y + radius * py, z + radius * pz
            if not any(j != i and (sx - atoms[j][2]) ** 2 + (sy - atoms[j][3]) ** 2 + (sz - atoms[j][4]) ** 2 < expanded[j] ** 2 for j in nearby):
                available += 1
        output[atom[0]] += (available / points) * 4 * math.pi * radius * radius
    return output


def alphafold_rsa(accession):
    records = json.loads(get(ALPHAFOLD + urllib.parse.quote(accession)))
    if not records:
        raise RuntimeError("No AlphaFold DB prediction was returned")
    model = records[0]
    pdb_url = model.get("pdbUrl")
    if not pdb_url:
        raise RuntimeError("AlphaFold DB response did not include a PDB URL")
    atoms = pdb_atoms(get(pdb_url, "text/plain"))
    sasa = residue_sasa(atoms)
    names = {atom[0]: atom[1] for atom in atoms}
    rsa = {number: 100 * area / MAX_ASA[names[number]] for number, area in sasa.items() if names[number] in MAX_ASA}
    return rsa, pdb_url


def exposure(value):
    if value is None:
        return "Unavailable"
    if value >= 35:
        return "Highly exposed"
    if value >= 20:
        return "Exposed"
    return "Buried / core"


def entry_name(entry):
    protein = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "-")
    genes = entry.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "-") if genes else "-"
    return protein, gene


def report(token, include_rsa=True):
    entry = entry_for(token)
    accession = entry["primaryAccession"]
    protein, gene = entry_name(entry)
    sequence = entry.get("sequence", {})
    rsa, pdb_url, rsa_error = {}, None, None
    if include_rsa:
        try:
            rsa, pdb_url = alphafold_rsa(accession)
        except RuntimeError as err:
            rsa_error = str(err)
    lines = [f"# Domain Architecture & Structural Accessibility for {accession} ({entry.get('uniProtkbId', '-')})",
             f"**Protein:** {protein} | **Gene:** {gene} | **Total length:** {sequence.get('length', '-')} aa",
             f"**Retrieved:** {dt.date.today().isoformat()} | **UniProt:** {UNIPROT}/{accession}.json"]
    if include_rsa:
        lines.append(f"**AlphaFold model:** {pdb_url or 'Unavailable'}")
        if rsa_error:
            lines.append(f"**RSA warning:** {rsa_error}")
    lines += ["", "## Defined domains, regions, and motifs",
              "| Type | Start | End | Length | Description | Structural RSA (%) | Exposure | Membrane context | Feature ID |",
              "| --- | ---: | ---: | ---: | --- | ---: | --- | --- | --- |"]
    for row in selected_features(entry):
        values = [rsa[p] for p in range(row["start"], row["end"] + 1) if p in rsa]
        average = sum(values) / len(values) if values else None
        membrane = "Membrane-associated (annotation)" if row["type"] in {"Transmembrane", "Intramembrane"} else "Not membrane"
        rsa_text = f"{average:.1f}" if average is not None else "-"
        lines.append(f"| {row['type']} | {row['start']} | {row['end']} | {row['end'] - row['start'] + 1} | {row['description']} | {rsa_text} | {exposure(average)} | {membrane} | {row['id']} |")
    lines += ["", "Structural RSA is calculated from an AlphaFold predicted structure (Shrake-Rupley; 1.4 Å probe; 240 points). Membrane context is a UniProt annotation and does not modify the calculated value."]
    return "\n".join(lines), entry


def input_tokens(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        dialect = csv.excel_tab if path.lower().endswith((".tsv", ".tab")) else csv.excel
        for index, row in enumerate(csv.reader(handle, dialect=dialect), start=1):
            candidates = [value.strip().upper() for value in row[:2] if value.strip()]
            token = next((value for value in candidates if ACCESSION.match(value) or "_" in value), None)
            if token:
                yield index, token


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-a", "--accession", help="UniProt accession or entry name")
    group.add_argument("-i", "--input", help="CSV, TSV, or text input (accession in column 1 or 2)")
    group.add_argument("-s", "--search", help="UniProtKB search query")
    parser.add_argument("-f", "--format", choices=("domains", "summary", "fasta"), default="domains")
    parser.add_argument("-o", "--output", help="Write report to this path instead of stdout")
    parser.add_argument("--no-rsa", action="store_true", help="Skip AlphaFold retrieval and structural RSA calculation")
    args = parser.parse_args()
    if args.search:
        query = urllib.parse.urlencode({"query": args.search, "format": "tsv", "fields": "accession,id,protein_name,gene_names,organism_name", "size": 25})
        result = get(f"{UNIPROT}/search?{query}", "text/tab-separated-values")
    elif args.input:
        chunks = []
        for line, token in input_tokens(args.input):
            try:
                text, _ = report(token, not args.no_rsa)
                chunks.append(text)
            except RuntimeError as err:
                chunks.append(f"# {token}\n**Input line {line} failed:** {err}")
        result = "\n\n---\n\n".join(chunks) or "No accession-like values found in the first two columns."
    else:
        entry = entry_for(args.accession)
        if args.format == "fasta":
            result = get(f"{UNIPROT}/{entry['primaryAccession']}.fasta", "text/plain")
        elif args.format == "summary":
            result = json.dumps(entry, indent=2)
        else:
            result, _ = report(entry["primaryAccession"], not args.no_rsa)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(result + ("" if result.endswith("\n") else "\n"))
    else:
        print(result)


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(2)
