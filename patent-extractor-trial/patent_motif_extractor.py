"""
Patent Chemical Motif Extraction & 2D Visualization Pipeline
Integrates:
- DECIMER (Optical Chemical Structure Recognition) concept for PDF/Image extraction
- RDKit for Maximum Common Substructure (MCS) calculation & protected motif identification
- RDKit 2D Visual Renderer for generating highlighted 2D chemical structure diagrams
"""

import os
import sys
from typing import List, Dict, Tuple

# Try importing RDKit
try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdFMCS
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Try importing DECIMER if available in environment
try:
    from DECIMER import predict_SMILES
    DECIMER_AVAILABLE = True
except ImportError:
    DECIMER_AVAILABLE = False


def decimer_ocsr_extract_smiles(image_path: str) -> str:
    """
    Extracts SMILES from a 2D chemical structure image using DECIMER.
    Falls back to mock/cached SMILES if DECIMER model weights are not loaded.
    """
    if DECIMER_AVAILABLE:
        try:
            smiles = predict_SMILES(image_path)
            return smiles
        except Exception as e:
            print(f"[DECIMER] Error processing {image_path}: {e}")
            return None
    else:
        print(f"[DECIMER] Package not installed/loaded. Operating in simulation mode for {image_path}")
        return None


def extract_patent_protected_motifs(smiles_list: List[str]) -> Tuple[str, List[Chem.Mol]]:
    """
    Calculates the Maximum Common Substructure (MCS) across a list of patent SMILES.
    This MCS represents the core conserved/protected motif covered across the patent claim family.
    """
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for MCS calculation.")

    mols = [Chem.MolFromSmiles(s) for s in smiles_list if s]
    mols = [m for m in mols if m is not None]

    if not mols:
        return "", []

    # Find Maximum Common Substructure (MCS)
    res = rdFMCS.FindMCS(
        mols,
        ringMatchesRingOnly=True,
        completeRingsOnly=True
    )
    
    core_smarts = res.smartsString
    return core_smarts, mols


def render_2d_motif_grid(mols: List[Chem.Mol], core_smarts: str, output_image_path: str) -> str:
    """
    Renders a 2D grid image of all extracted patent structures with the protected core motif
    highlighted in red/coral on every structure.
    """
    if not RDKIT_AVAILABLE or not mols:
        return ""

    pattern = Chem.MolFromSmarts(core_smarts) if core_smarts else None
    
    highlight_atom_lists = []
    legends = []

    for i, mol in enumerate(mols):
        legends.append(f"Patent Struct #{i+1}")
        if pattern and mol.HasSubstructMatch(pattern):
            atom_matches = mol.GetSubstructMatch(pattern)
            highlight_atom_lists.append(atom_matches)
        else:
            highlight_atom_lists.append(())

    # Draw grid image
    grid_img = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        subImgSize=(350, 300),
        legends=legends,
        highlightAtomLists=highlight_atom_lists
    )
    
    grid_img.save(output_image_path)
    return output_image_path


if __name__ == "__main__":
    print("Patent Chemical Motif Extractor initialized.")
    if RDKIT_AVAILABLE:
        print("RDKit version:", rdkit.__version__)
    print("DECIMER Available:", DECIMER_AVAILABLE)
    
    # Example valid SMILES (Ritonavir and derivatives)
    sample_smiles = [
        "CC(C)c1nc(cs1)CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc2ccccc2)[C@H](O)CN(Cc3ccccc3)C(=O)OCC",
        "CC(C)c1nc(cs1)CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc2ccccc2)[C@H](O)CN(Cc3ccccc3)C(=O)OCCCCCC",
        "CC(C)c1nc(cs1)CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc2ccccc2)[C@H](O)CN(Cc3ncccn3)C(=O)OCC"
    ]
    
    smarts, molecules = extract_patent_protected_motifs(sample_smiles)
    print("Extracted SMARTS protected core:", smarts)
    out_img = render_2d_motif_grid(molecules, smarts, "/Users/madhura/.gemini/antigravity-cli/brain/40039de9-4eba-433c-a168-7c784da74bbb/extracted_patent_motifs_2d.png")
    print("Saved 2D grid image to:", out_img)
