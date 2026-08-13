# Apixaban Patent Exclusivity & Motif Analysis Report

This report summarizes the computational analysis of the patent-protected structural motifs of **Apixaban (Eliquis)** and outlines the Loss of Exclusivity (LOE) timeline and design-around feasibility.

---

## 1. Drug Profile & Exclusivity Timeline

- **Brand Name**: Eliquis®
- **Active Ingredient**: Apixaban
- **Drug Class / Mechanism**: Factor Xa Inhibitor (Anticoagulant)
- **Primary Patent**: US 6,967,208 (BMS / Pfizer)
- **Calculated Primary Expiration**: November 2026 (including Patent Term Extension).
- **Pediatric Exclusivity Extension**: Expiring May 2028.
- **Formulation Patent**: US 9,326,945 (Expires November 2031).

---

## 2. Conserved Patent-Protected Core Motif

The Maximum Common Substructure (MCS) algorithm identified the conserved, legally protected core motif shared across Apixaban and its claimed synthetic analogs:

- **Conserved SMARTS Scaffold**:
  `[#6&!R]-&!@[#8&!R]-&!@[#6]1:&@[#6]:&@[#6]:&@[#6](:&@[#6]:&@[#6]:&@1)-&!@[#7]1:&@[#6]2:&@[#6](:&@[#6](:&@[#7]:&@1)-&!@[#6&!R](=&!@[#8&!R])-&!@[#7&!R])-&@[#6]-&@[#6]-&@[#7](-&@[#6]-&@2=&!@[#8&!R])-&!@[#6]1:&@[#6]:&@[#6]:&@[#6](:&@[#6]:&@[#6]:&@1)-&!@[#7&R]`

This represents:
1. **Central Dihydropyrazolo-pyridinone Core**: The central bicyclic system essential for anchoring inside the Factor Xa active site.
2. **Methoxyphenyl Group**: Fits into the S1 pocket of the enzyme.
3. **Oxopiperidin-phenyl Ring Linker**: Anchors the terminal group in the S4 pocket.

---

## 3. 2D Visual Motif Grid

The 2D visual grid below highlights the shared, patent-protected core scaffold in red, leaving the variable points (such as the terminal amide or modified lactam rings) uncolored:

![Apixaban Protected Motif 2D Grid](file:///Users/madhura/trial-skills/patent_extraction/apixaban_motifs_2d.png)

---

## 4. Opportunity Summary & Design-Around Feasibility

### A. Generic Small Molecule Entry
- **Primary Compound Patent (US 6,967,208)**: Expiring late 2026. A generic entry is highly viable post-2026/2028 depending on target jurisdictions.
- **Formulation Patent Challenge (US 9,326,945)**: BMS/Pfizer successfully defended this patent in recent litigations, delaying generic entry in the US to 2028 under settlement agreements. Generic developers must use non-infringing formulation strategies (e.g. avoiding the claimed solid dispersion dissolution profiles) or launch post-settlement windows.

### B. Structural Modification (NCE Derivative / Bioisostere)
- **Modifiable Positions**:
  - **Terminal Amide**: The carboxamide group can be converted to esters or alternative substituted amides (Variant A).
  - **Lactam Ring**: The 2-piperidone ring can be substituted with oxomorpholine (Variant B) or other saturated heterocycles without destroying Factor Xa affinity, offering potential pathways to novel, non-infringing analogs.
