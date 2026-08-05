# PRSS55 intact-protein substrate analysis

This project prioritizes predicted 8-residue PRSS55 substrate windows for
wet-lab testing. It is a deterministic structural-analysis workflow: it does
**not** train, modify, or recalibrate the upstream cleavage-prediction model.

The workflow validates each predicted sequence against UniProt, applies
precursor/topology exclusions, calculates AlphaFold-model-derived structural
features, assigns three-state secondary structure, and produces an auditable
shortlist ordered by the supplied `prob_cleavage` value.

## Scientific guardrails

- A `prob_cleavage` value ranks the upstream predictions; it is not evidence
  that an intact human protein is cleaved.
- The workflow assigns a site `ELIGIBLE`, `REVIEW`, or `EXCLUDE`; it never calls
  a confirmed substrate.
- RSA and secondary structure are derived from one AlphaFold model. They do not
  account for native complexes, glycans, membranes, processing heterogeneity, or
  conformational ensembles.
- Verify that the predicted P1/P1′ cleavage convention is known before using
  bond-level features. In the current PRSS55 dataset, the central arginine is
  treated as P1 (`Arg↓X`), therefore `--p1-offset 4` is used.

## Install

Python 3.13 or later is expected. The structural and secondary-structure
dependencies are installed project-locally:

```bash
pip3 install --target .vendor -r requirements.txt
```

The first extraction run sends each UniProt accession to the public UniProt and
AlphaFold endpoints to retrieve annotations and structures. UniProt JSON is
cached under `cache/uniprot/`; derived scores and result tables are computed and
stored locally.

## Input

The raw prediction table must contain these columns:

```text
seq_id  protein_name  window_idx  start_pos  end_pos  8mer_sequence
prob_cleavage  prediction
```

`seq_id` must be in a UniProt-style format such as `sp|Q6PDA7|SG11A_HUMAN`.
Coordinates must be 1-based and refer to the same canonical UniProt sequence
identified by the accession. See `data/prss55_predictions.raw.tsv` for an
example.

## Run the workflow

Run from the project root.

```bash
# 1. Normalize the model output and record the known P1 position.
python3 scripts/prepare_prss55_predictions.py \
  --input data/prss55_predictions.raw.tsv \
  --output data/prss55_predicted_sites.csv \
  --p1-offset 4

# 2. Validate sequence coordinates, retrieve annotations, compute RSA/pLDDT,
#    and assign DSSP-compatible secondary structure.
PYTHONPATH=.vendor python3 scripts/extract_prss55_site_features.py \
  --input data/prss55_predicted_sites.csv \
  --output results/prss55_site_features.tsv \
  --with-rsa

# 3. Apply deterministic exclusions and sort the remaining rows by prob_cleavage.
python3 scripts/triage_prss55_sites.py \
  --features results/prss55_site_features.tsv \
  --output results/prss55_triaged_sites.tsv
```

## Outputs

`results/prss55_site_features.tsv` contains one row per 8-mer with:

- sequence-mapping status and mature-protein context (`Signal`, `Propeptide`,
  `Peptide`, `Chain`);
- `domain_annotations`: overlap with UniProt `Domain`, `Region`, `Repeat`,
  `Motif`, `Zinc finger`, or `Coiled coil` features, retaining type and
  description;
- `family_annotations`: InterPro/Pfam membership, reported separately because
  those records may not provide residue boundaries;
- P1/P1′ and window RSA, plus mean pLDDT;
- PyDSSP three-state annotation: `H` alpha helix, `E` beta strand, `C` loop or
  other; and fractions of each state over the window;
- PTM/glycosylation, disulfide, topology, and evidence-gap flags;
- exact AlphaFold model URL used for the calculation.

`results/prss55_triaged_sites.tsv` adds:

- `ELIGIBLE`: no automatic sequence, processing, or membrane conflict;
- `REVIEW`: structurally possible but requiring human interpretation;
- `EXCLUDE`: sequence mismatch, signal/propeptide overlap, or membrane overlap.

## Required review before selecting wet-lab targets

1. Confirm the P1/P1′ convention in the prediction labels.
2. Inspect every `EXCLUDE` and `REVIEW` row; do not remove a processing flag
   solely because the model score is high.
3. Check that each retained site is physically reachable in the specific human
   sample type and relevant PRSS55-accessible compartment.
4. Treat disulfide overlap, nearby glycosylation, low AlphaFold confidence, and
   missing topology as reasons for manual review.
5. Archive the raw input, generated TSVs, UniProt cache, and database/model URLs
   with every experimental batch.

Further method details are in
[`docs/prss55_intact_protein_substrate_prioritization.md`](docs/prss55_intact_protein_substrate_prioritization.md)
and [`docs/prss55_pipeline_usage.md`](docs/prss55_pipeline_usage.md).
