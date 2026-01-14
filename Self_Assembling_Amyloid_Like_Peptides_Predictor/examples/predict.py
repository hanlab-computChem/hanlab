"""
Peptide prediction script using UniMol model.

This script predicts amyloid-like peptide properties (AP and SHB scores)
for given peptide sequences using a pre-trained UniMol model.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from unimol_tools import MolPredict


def get_project_root() -> Path:
    """Get the project root directory (parent of examples directory)."""
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent


def check_confgen_available() -> bool:
    """Check if confgen command is available in PATH.
    
    Returns:
        True if confgen is available, False otherwise.
    """
    try:
        subprocess.run(
            ["confgen", "--version"],
            check=True,
            capture_output=True,
            text=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def process_mol_with_rdkit(smiles: str, peptide_name: str) -> List[Tuple[List[str], Any]]:
    """Process a single molecule using RDKit to generate 3D conformers.
    
    Args:
        smiles: SMILES string of the molecule.
        peptide_name: Name/sequence of the peptide.
        
    Returns:
        List of tuples containing (atoms, coordinates) for each conformer.
        
    Raises:
        RuntimeError: If molecule cannot be created or processed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"Failed to create molecule from SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)
    params = AllChem.EmbedParameters()
    params.pruneRmsThresh = 2.0
    
    # Generate multiple conformers
    num_confs = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
    if num_confs == 0:
        raise RuntimeError(f"Failed to generate conformers for {peptide_name}")
    
    # Optimize conformers using MMFF
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise RuntimeError(f"Failed to get MMFF properties for {peptide_name}")
    
    for conf_id in range(mol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    
    # Calculate energies and sort
    energies = [
        (conf.GetId(), AllChem.MMFFGetMoleculeForceField(
            mol, props, confId=conf.GetId()
        ).CalcEnergy())
        for conf in mol.GetConformers()
    ]
    energies.sort(key=lambda x: x[1])
    
    # Select top 10 conformers
    selected_conf_ids = [conf_id for conf_id, _ in energies[:10]]
    
    # Extract atoms and coordinates
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    results = []
    for conf_id in selected_conf_ids:
        conf = mol.GetConformer(conf_id)
        coords = conf.GetPositions()
        results.append((atoms, coords))
    
    print(f"{peptide_name} processed with RDKit")
    return results


def generate_3d_coordinates_with_rdkit(
    peptides: List[str],
    seq2smi: Dict[str, str]
) -> Dict[str, List]:
    """Generate 3D coordinates for peptides using RDKit.
    
    Args:
        peptides: List of peptide sequences to process.
        seq2smi: Dictionary mapping sequences to SMILES strings.
        
    Returns:
        Dictionary containing atoms, sequences, and coordinates.
        
    Raises:
        RuntimeError: If no valid molecules are generated.
    """
    data_to_pred = {"atoms": [], "seqs": [], "coordinates": []}
    
    for pep in peptides:
        pep_upper = pep.upper()
        if pep_upper not in seq2smi:
            print(f"Warning: Peptide [{pep_upper}] is not supported yet. Skipping.")
            continue
        
        smiles = seq2smi[pep_upper]
        try:
            conformers = process_mol_with_rdkit(smiles, pep_upper)
            for atoms, coords in conformers:
                data_to_pred["atoms"].append(atoms)
                data_to_pred["seqs"].append(pep_upper)
                data_to_pred["coordinates"].append(coords)
        except RuntimeError as e:
            print(f"Warning: Failed to process {pep_upper}: {e}. Skipping.")
            continue
    
    if not data_to_pred["seqs"]:
        raise RuntimeError("No valid molecules were generated using RDKit.")
    
    return data_to_pred


def load_sequence_to_smiles_mapping(seq2smi_path: Path) -> Dict[str, str]:
    """Load sequence to SMILES mapping from JSON file.
    
    Args:
        seq2smi_path: Path to the JSON file containing sequence to SMILES mapping.
        
    Returns:
        Dictionary mapping peptide sequences to SMILES strings.
        
    Raises:
        FileNotFoundError: If the mapping file does not exist.
        json.JSONDecodeError: If the file is not valid JSON.
    """
    if not seq2smi_path.exists():
        raise FileNotFoundError(f"Sequence to SMILES mapping file not found: {seq2smi_path}")
    
    with open(seq2smi_path, "r", encoding="utf-8") as f:
        return json.load(f)


def generate_3d_coordinates(
    peptides: List[str],
    seq2smi: Dict[str, str],
    work_dir: Path
) -> Dict[str, List]:
    """Generate 3D coordinates for peptides using confgen or RDKit as fallback.
    
    Args:
        peptides: List of peptide sequences to process.
        seq2smi: Dictionary mapping sequences to SMILES strings.
        work_dir: Working directory for temporary files.
        
    Returns:
        Dictionary containing atoms, sequences, and coordinates.
        
    Raises:
        RuntimeError: If confgen command fails or no valid molecules are generated.
    """
    # Check if confgen is available
    if check_confgen_available():
        print("Using confgen to generate 3D coordinates...")
        return _generate_3d_coordinates_with_confgen(peptides, seq2smi, work_dir)
    else:
        print("confgen not found. Using RDKit as fallback...")
        return generate_3d_coordinates_with_rdkit(peptides, seq2smi)


def _generate_3d_coordinates_with_confgen(
    peptides: List[str],
    seq2smi: Dict[str, str],
    work_dir: Path
) -> Dict[str, List]:
    """Generate 3D coordinates for peptides using confgen.
    
    Args:
        peptides: List of peptide sequences to process.
        seq2smi: Dictionary mapping sequences to SMILES strings.
        work_dir: Working directory for temporary files.
        
    Returns:
        Dictionary containing atoms, sequences, and coordinates.
        
    Raises:
        RuntimeError: If confgen command fails or no valid molecules are generated.
    """
    # Prepare SMILES file
    smi_file = work_dir / "pep_to_pred.smi"
    sdf_file = work_dir / "pep_to_pred.sdf"
    
    smilist = []
    for pep in peptides:
        pep_upper = pep.upper()
        if pep_upper not in seq2smi:
            print(f"Warning: Peptide [{pep_upper}] is not supported yet. Skipping.")
            continue
        smilist.append(f"{seq2smi[pep_upper]}\t{pep_upper}")
    
    if not smilist:
        raise ValueError("No valid peptides found to process.")
    
    # Write SMILES file
    with open(smi_file, "w", encoding="utf-8") as f:
        f.write("\n".join(smilist))
    
    # Run confgen to generate 3D coordinates
    confgen_cmd = [
        "confgen",
        "-i", str(smi_file),
        "-o", str(sdf_file),
        "-C", "LARGE_SET_DIVERSE",
        "-v", "ERROR",
        "-T", "60",
        "-n", "10",
        "-r", "2.0"
    ]
    
    try:
        result = subprocess.run(
            confgen_cmd,
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"confgen command failed: {e.stderr}") from e
    except FileNotFoundError:
        raise RuntimeError("confgen command not found. Please ensure it is installed and in PATH.")
    
    # Parse SDF file
    data_to_pred = {"atoms": [], "seqs": [], "coordinates": []}
    mol_supplier = Chem.SDMolSupplier(str(sdf_file))
    
    for mol in mol_supplier:
        if mol is None:
            continue
        data_to_pred["atoms"].append([atom.GetSymbol() for atom in mol.GetAtoms()])
        data_to_pred["seqs"].append(mol.GetProp("_Name"))
        data_to_pred["coordinates"].append(mol.GetConformer().GetPositions())
    
    if not data_to_pred["seqs"]:
        raise RuntimeError("No valid molecules were generated from the SDF file.")
    
    return data_to_pred


def predict_with_unimol(
    data_to_pred: Dict[str, List],
    model_path: Path
) -> Dict[str, Tuple[float, float]]:
    """Run UniMol prediction on the prepared data.
    
    Args:
        data_to_pred: Dictionary containing atoms, sequences, and coordinates.
        model_path: Path to the trained model directory.
        
    Returns:
        Dictionary mapping peptide sequences to predicted scores (AP, SHB).
        
    Raises:
        FileNotFoundError: If the model directory does not exist.
    """
    if not model_path.exists():
        raise FileNotFoundError(f"Model directory not found: {model_path}")
    
    clf = MolPredict(load_model=str(model_path))
    predictions = clf.predict(data=data_to_pred).tolist()
    
    # Aggregate predictions for peptides with multiple conformers
    temp_result: Dict[str, List[List[float]]] = {}
    for pep, score in zip(data_to_pred["seqs"], predictions):
        if pep not in temp_result:
            temp_result[pep] = []
        temp_result[pep].append(score)
    
    # Calculate average scores for each peptide
    result_score: Dict[str, Tuple[float, float]] = {}
    for pep, scores in temp_result.items():
        avg_ap = sum(s[0] for s in scores) / len(scores)
        avg_shb = sum(s[1] for s in scores) / len(scores)
        result_score[pep] = (avg_ap, avg_shb)
    
    return result_score


def save_results(result_score: Dict[str, Tuple[float, float]], output_path: Path) -> None:
    """Save prediction results to CSV file.
    
    Args:
        result_score: Dictionary mapping peptides to predicted scores.
        output_path: Path to the output CSV file.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Peptide,Pred_AP,Pred_SHB\n")
        for pep, (ap_score, shb_score) in sorted(result_score.items()):
            f.write(f"{pep},{ap_score:.4f},{shb_score:.4f}\n")
    
    print(f"Results saved to: {output_path}")


def main():
    """Main function to run peptide prediction pipeline."""
    parser = argparse.ArgumentParser(
        description="Predict amyloid-like peptide properties using UniMol model"
    )
    parser.add_argument(
        "-p", "--peptides",
        type=str,
        required=True,
        nargs="+",
        help="List of peptide sequences to predict (space-separated)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output CSV file path (default: predict_result.csv in current directory)"
    )
    args = parser.parse_args()
    
    # Setup paths
    project_root = get_project_root()
    data_dir = project_root / "data"
    model_dir = project_root / "trained_model"
    work_dir = Path(__file__).resolve().parent
    
    seq2smi_path = data_dir / "seq2smi_1to5.json"
    output_path = Path(args.output) if args.output else work_dir / "predict_result.csv"
    
    try:
        # Load sequence to SMILES mapping
        print("Loading sequence to SMILES mapping...")
        seq2smi = load_sequence_to_smiles_mapping(seq2smi_path)
        
        # Normalize peptide sequences
        peptides = [pep.upper() for pep in args.peptides]
        
        # Generate 3D coordinates
        print("Generating 3D coordinates...")
        data_to_pred = generate_3d_coordinates(peptides, seq2smi, work_dir)
        
        # Run prediction
        print("Running UniMol prediction...")
        result_score = predict_with_unimol(data_to_pred, model_dir)
        
        # Save results
        save_results(result_score, output_path)
        
        print(f"Successfully predicted {len(result_score)} peptides.")
        
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()