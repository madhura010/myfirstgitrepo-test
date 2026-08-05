# Computational prioritization of intact-protein substrates for human PRSS55

## Purpose

Rank predicted PRSS55 8-mer sites in human proteins for follow-up testing. This
is a **cleavability-prioritization** workflow, not a claim that a protein has
been cleaved. Computational analysis cannot demonstrate a proteolytic
neo-terminus or an exact cleavage event.

## Scope and assumptions

- An upstream prediction tool supplies predicted 8-mers, their protein
  coordinates, and a `prob_cleavage` score. For the present PRSS55 dataset, the
  central arginine (position 4) is P1 and the proposed bond is Arg↓P1′.
- The objective is cleavage of intact proteins, not merely hydrolysis of short
  peptides.
- Co-expression in the relevant reproductive setting has already been applied
  as an eligibility filter and is not re-scored here.
- Each candidate is analysed as a *site* (P4–P4′), not as a protein-level
  property.

## Biological constraints specific to PRSS55

Human PRSS55 is testis-enriched and annotated as a probable S1 serine protease.
Its membrane topology is uncertain; the human entry notes possible GPI anchoring
or a type-I membrane topology. Mouse data support a GPI-anchored protein on the
sperm/acrosomal surface. Consequently, a candidate must be present on the side
of a membrane or in the compartment that PRSS55 can physically access. Do not
assume that every co-expressed intracellular protein is eligible.

Sources: [UniProt PRSS55 (Q6UWB4)](https://www.uniprot.org/uniprotkb/Q6UWB4-2);
[Shang et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6208766/).

## Workflow

### 1. Preserve and interpret the prediction input correctly

Before applying structural ranking, verify that each window maps to the correct
protein sequence and preserve the score direction and range.

- Store each hit as an 8-residue window with its start/end coordinates.
- Use `prob_cleavage` only to sort candidates after structural exclusions; it
  is not a probability that the intact protein is cleaved.
- Do not invent a P1-P1′ bond for future datasets. An 8-mer has seven possible
  peptide bonds. Here, use the stated central-arginine convention (P1 offset 4)
  and document it in every result.

### 2. Map every hit to the correct protein form

For each hit, record:

- UniProt accession, reviewed isoform, and sequence version;
- 8-mer sequence and residue coordinates;
- P1–P1′ coordinate and flanking sequence;
- whether the site is removed by a signal peptide, pro-peptide, or other known
  processing event;
- whether the same site occurs in multiple isoforms.

Exclude hits that cannot be mapped unambiguously or do not occur in the mature
protein form being evaluated.

### 3. Apply hard accessibility and topology exclusions

Discard or flag candidates when the scissile bond:

- lies in a transmembrane segment or signal peptide;
- is on the membrane side inaccessible to PRSS55;
- is in a cytosolic/nuclear region incompatible with the intended PRSS55
  compartment;
- is buried in the core of a folded domain;
- has no usable structure or disorder evidence and cannot be assessed.

Use UniProt topology/processing annotations first. Use a high-quality
experimental structure when available; otherwise use an AlphaFold model, marked
as predicted structural evidence.

### 4. Compute site-level structural features

Calculate all structural features over the 8-mer window and at the defined
P1–P1′ bond (Arg position 4 followed by position 5).

| Feature | Calculation | Interpretation |
| --- | --- | --- |
| Window accessibility | Mean and minimum SASA/RSA across the 8-mer | Site-level evidence of exposure; do not use protein-wide or domain-average SASA alone. |
| Bond-specific accessibility | Per-residue SASA/RSA at P1 and P1′, once the bond is defined | Do not compute or infer this from the centre of the 8-mer. |
| Structural depth/contact density | Local atom/residue contacts or depth below protein surface | Penalize grooves and cores that are nominally solvent exposed but sterically constrained. |
| Secondary structure | Secondary structure at least across P2–P2′ | Favor coil, turn, or accessible loop; penalize stable alpha-helix and beta-strand. |
| Flexibility | Disorder prediction plus AlphaFold pLDDT where applicable | Supports, but does not prove, transient access. Low pLDDT is not direct evidence of cleavage. |
| Domain context | Domain/core, disulfide-rich region, catalytic site, or binding region | Penalize sites where cleavage is structurally implausible or disruptive without evidence. |
| Interface occlusion | Experimental complex structures or credible complex/interface annotations | Penalize monomer-exposed sites that may be buried in the native assembly. |
| PTM/glycan shielding | Known glycosylation/PTM annotations and nearby N-X-S/T sequons | Flag sites likely shielded in the native protein. |
| Structural confidence | Structure source, model confidence, coverage | Record uncertainty explicitly rather than forcing a precise score. |

### 5. Assess exposure robustness

Avoid treating one static AlphaFold structure as ground truth. Assign an exposure
robustness label for each site:

- **High:** exposed in an experimental structure or consistently exposed across
  multiple relevant structural observations.
- **Medium:** exposed in one confident predicted structural region.
- **Low:** model is low confidence, site may be an interface, or accessibility
  differs across plausible structures.

Whole-protein PRSS55–substrate docking is not a primary filter: unconstrained
docking across many targets tends to create precise-looking but unvalidated
poses. It may be used only for qualitative inspection of the final few sites.

### 6. Rank candidates transparently

Use hard exclusions before any score. For surviving sites, rank with a
documented composite:

```text
priority =
    hard exclusions (sequence mismatch, processing overlap, membrane conflict)
  → manual-review flags (topology, glycan/PTM, disulfide, missing structure)
  → sort remaining sites by supplied prob_cleavage
```

Until positive and negative experimental labels exist, this is a prioritization
workflow, not a calibrated probability. Use statuses:

- **ELIGIBLE:** no automatic sequence, processing, or membrane conflict.
- **REVIEW:** no hard conflict, but topology/structure/PTM evidence needs human
  interpretation.
- **EXCLUDE:** incompatible sequence mapping, processing annotation, or
  transmembrane/intramembrane overlap.

### 7. Produce a reproducible candidate table

One row per predicted site:

```text
accession | isoform | protein name | 8-mer | P1-P1′ coordinate |
prob_cleavage | P1/P1′ SASA | min/mean window SASA |
secondary structure | disorder/pLDDT | topology status |
PTM/interface flags | structure source/confidence |
exposure robustness | composite score | tier | rationale
```

Keep the raw feature values, tool/database versions, structure identifiers, and
exclusion rationale. This permits later calibration using experimental outcomes.

## What experimental testing will establish later

Testing the top-ranked sites should determine whether cleavage occurs in the
intact substrate and, if so, which bond is cut. Until then, the pipeline reports
ranked hypotheses, not confirmed PRSS55 substrates.
