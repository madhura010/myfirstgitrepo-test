# PRSS55 site-prioritization tools

## What the tools do

`scripts/extract_prss55_site_features.py` validates each predicted 8-mer against
the current UniProt sequence, retrieves maturity/topology/PTM annotations, and
optionally calculates per-residue RSA from an AlphaFold model. It does **not**
infer secondary structure, protein-complex interfaces, or a cleavage probability
without separately supplied evidence.

The output distinguishes a UniProt `Chain` (the retained mature protein) from a
`Signal` or `Propeptide` (a region removed during processing), and separately
records an annotated mature `Peptide`. This distinction is essential for small
secreted proteins such as defensins.

`domain_annotations` contains residue-bounded UniProt Domain, Region, Repeat,
Motif, Zinc finger, and Coiled coil features. Each entry retains its feature
type plus description. `family_annotations` contains InterPro/Pfam memberships,
which are useful context but are not reported as site coordinates unless the
source actually supplies boundaries.

With the project-local PyDSSP dependency installed, the output also includes
`ss_window_pattern` and `ss_p1`/`ss_p1prime`: `H` = alpha helix, `E` = beta
strand, and `C` = loop/other. The `ss_*_fraction` columns are the fraction of
the eight residues assigned to each state. They are assignments for the chosen
AlphaFold model, not an experimental ensemble measurement.

`scripts/triage_prss55_sites.py` is a deterministic rule-based analysis step.
It excludes sequence/processing/transmembrane conflicts, flags unresolved
evidence, then sorts eligible candidates by the supplied `prob_cleavage`. It
does not fit, retrain, calibrate, or create an additional ML model.

## Input schema

The supplied table can be normalized directly:

```bash
python3 scripts/prepare_prss55_predictions.py \
  --input data/prss55_predictions.raw.tsv \
  --output data/prss55_predicted_sites.csv \
  --p1-offset 4
```

For the current PRSS55 predictions, the pivotal central arginine is position 4,
so the proposed cleavage is `P4-P3-P2-Arg↓P1′-P2′-P3′-P4′`. The `--p1-offset 4`
argument records this explicit biological convention rather than guessing it.

For other input sources, start from `data/prss55_predicted_sites.template.csv`.
Required fields are:

| Field | Meaning |
| --- | --- |
| `accession` | UniProt accession for the candidate protein |
| `site_start` | 1-based start position of the 8-mer in that accession's sequence |
| `eight_mer` | Eight amino-acid one-letter codes |
| `p1_offset` | Optional 1-based P1 position within the 8-mer; leave blank when the exact bond is unknown |
| `model_score` | Supplied `prob_cleavage`; used only to sort structurally eligible candidates |

`source_isoform` and `source_id` are optional but strongly recommended.

## Run feature extraction

Install the documented optional structural dependency once when secondary-
structure assignment is required:

```bash
pip3 install --target .vendor -r requirements.txt
```

```bash
python3 scripts/extract_prss55_site_features.py \
  --input data/prss55_predicted_sites.csv \
  --output results/prss55_site_features.tsv \
  --with-rsa
```

The first run retrieves UniProt records into `cache/uniprot/`. Inspect every row
with a sequence mismatch, processing overlap, unknown topology, missing
structure, or incomplete structural coverage before proceeding.

## Run deterministic triage

```bash
python3 scripts/triage_prss55_sites.py \
  --features results/prss55_site_features.tsv \
  --output results/prss55_triaged_sites.tsv
```

`ELIGIBLE` means no automatic structural conflict was found; `REVIEW` means
manual structural/topology interpretation is required; `EXCLUDE` means a
sequence-processing or membrane conflict was found. None means confirmed
substrate.

## Required human review before ordering proteins

- Confirm P1–P1′ alignment and model-score direction/scale.
- Resolve mature-protein processing and membrane orientation for every retained
  candidate.
- Add experimental structure/complex/interface evidence where available.
- Review glycosylation, other PTMs, and disulfide annotations in the actual
  biological form. An annotated mature `Peptide` is supporting context, not an
  automatic exclusion.
- Record database release dates, model version, input checksum, and all manual
  inclusion/exclusion decisions.
