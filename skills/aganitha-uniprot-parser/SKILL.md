---
name: aganitha-uniprot-parser
description: Retrieve and analyze UniProt protein entries, feature/domain architecture, AlphaFold structures, and structural SASA/RSA. Use this skill whenever a user needs UniProt accessions or entry data combined with domain/feature coordinates, AlphaFold-model solvent exposure, SASA, RSA, or a batch accession report. Do not use it for generic FASTA manipulation or unrelated sequence analysis.
---

# UniProt Domain Architecture & Solvent Accessibility Skill

This skill provides automated tools to fetch, parse, and analyze protein domain architecture, feature annotations, and **3D Solvent Accessibility (RSA % / SASA)** directly from the **UniProt REST API** and **AlphaFold 3D structure models**.

## Core Features

- **Domain Architecture Parsing**: Extracts exact start and end coordinates, residue length, description, and feature ID for all domain-related feature types.
- **3D Solvent Accessibility Calculation (RSA %)**:
  - Automatically fetches predicted 3D atomic coordinates from **AlphaFold DB**.
  - Runs a standard **Shrake-Rupley algorithm** (1.4 Å probe; 240 sphere points) on the AlphaFold PDB coordinates to compute per-residue SASA and RSA %.
  - Separately reports a **membrane-context interpretation** for UniProt `Transmembrane` and `Intramembrane` features. This never overwrites the calculated structural RSA: a PDB-only calculation does not model a lipid bilayer.
  - Classifies calculated structural exposure as Highly exposed (≥35%), Exposed (20–<35%), or Buried/core (<20%). These are heuristic labels, not experimental measurements.
- **Broad Domain Feature Coverage**:
  - `Domain`, `Region`, `Repeat`, `Motif`, `Zinc finger`, `Coiled coil`, `Topological domain`, `Transmembrane`, `Intramembrane`
- **Batch Processing from File**: Reads text, TSV, or CSV files containing UniProt accession numbers in either the **1st or 2nd column** (automatically skipping headers).

## Bundled Helper Script Usage

The script is located at `scripts/fetch_uniprot.py`.

### 1. Single Accession Domain & Solvent Accessibility Parsing

```bash
# Calculate Domain Architecture + 3D Solvent Accessibility (RSA %)
python3 <skill-dir>/scripts/fetch_uniprot.py -a P00533 -f domains

# Fast mode without 3D RSA calculation
python3 <skill-dir>/scripts/fetch_uniprot.py -a P00533 -f domains --no-rsa
```

### 2. Batch Domain & Accessibility Extraction from File

```bash
# Read accessions from the first or second column and output a report
python3 <skill-dir>/scripts/fetch_uniprot.py -i input_accessions.tsv -o domain_accessibility_report.md
```

### 3. Additional Formats & Search

```bash
# Search UniProtKB by gene or keyword
python3 <skill-dir>/scripts/fetch_uniprot.py -s "gene:TP53 AND organism_id:9606"

# Get FASTA sequence format
python3 <skill-dir>/scripts/fetch_uniprot.py -a P01308 -f fasta

# Get full entry summary
python3 <skill-dir>/scripts/fetch_uniprot.py -a P68871 -f summary
```

## Expected Domain & Solvent Accessibility Output

```markdown
# Domain Architecture & Solvent Accessibility for P00533 (EGFR_HUMAN)
**Protein:** Epidermal growth factor receptor | **Gene:** EGFR | **Total Length:** 1210 aa

## Defined Domains, Regions & Motifs
| Type | Start | End | Length | Description | Structural RSA (%) | Exposure | Membrane context | Feature ID |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Topological domain | 25 | 645 | 621 aa | Extracellular | **28.0%** | Exposed | Not membrane | - |
| Repeat | 75 | 300 | 226 aa | Approximate | **24.1%** | Exposed | Not membrane | - |
| Transmembrane | 646 | 668 | 23 aa | Helical | **12.0%** | Buried / Core | Membrane-associated | - |
| Topological domain | 669 | 1210 | 542 aa | Cytoplasmic | **44.2%** | Highly Exposed | Not membrane | - |
| Domain | 712 | 979 | 268 aa | Protein kinase | **23.1%** | Exposed | Not membrane | - |
| Region | 1097 | 1137 | 41 aa | Disordered | **74.3%** | Highly Exposed | Not membrane | - |
```

## Workflow Instructions for Agents

1. When the user asks for domain details or solvent accessibility for a protein (e.g. `P00533`, `P04637`, `HBB_HUMAN`):
   - Run `fetch_uniprot.py -a <accession> -f domains` to fetch AlphaFold 3D model and output domain RSA percentages.
2. When the user passes a file containing UniProt accessions:
   - Use `fetch_uniprot.py -i <file_path>` to generate a batch solvent accessibility report.

## Method and limitations

- The report records the UniProt accession, AlphaFold model URL, and retrieval date. Treat calculated values as model-derived estimates.
- RSA is SASA normalized with the Tien *et al.* maximum-ASA reference values included in the script. Residues without a reference value are omitted from the RSA average.
- An AlphaFold model can be unavailable or can have incomplete residue coverage. The script preserves feature rows and marks the RSA value unavailable instead of inventing a value.
- Batch input is intentionally conservative: it checks columns 1 and 2 for valid-looking accessions or UniProt entry names, skips obvious headers, and reports rows that cannot be resolved.
