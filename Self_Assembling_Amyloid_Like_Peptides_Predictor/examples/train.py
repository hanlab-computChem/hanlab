"""
Training script for UniMol model on peptide data.

This script trains a UniMol model for multilabel regression task
to predict amyloid-like peptide properties.
"""

import argparse
import pickle
import sys
from pathlib import Path
from typing import Any

from unimol_tools import MolTrain


def get_project_root() -> Path:
    """Get the project root directory (parent of examples directory)."""
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent


def load_training_data(data_path: Path) -> Any:
    """Load training data from pickle file.
    
    Args:
        data_path: Path to the pickle file containing training data.
        
    Returns:
        Loaded training data.
        
    Raises:
        FileNotFoundError: If the data file does not exist.
        pickle.UnpicklingError: If the file cannot be unpickled.
    """
    if not data_path.exists():
        raise FileNotFoundError(f"Training data file not found: {data_path}")
    
    try:
        with open(data_path, "rb") as f:
            return pickle.load(f)
    except pickle.UnpicklingError as e:
        raise pickle.UnpicklingError(f"Failed to unpickle data file: {e}") from e


def train_model(
    data: Any,
    output_dir: Path,
    epochs: int = 500,
    batch_size: int = 128,
    learning_rate: float = 1e-4,
    early_stopping: int = 20,
    kfold: int = 4,
    gpu_id: int = 0,
    label_weight: list = None
) -> None:
    """Train UniMol model on peptide data.
    
    Args:
        data: Training data.
        output_dir: Directory to save trained model and outputs.
        epochs: Number of training epochs.
        batch_size: Batch size for training.
        learning_rate: Learning rate for optimizer.
        early_stopping: Number of epochs to wait before early stopping.
        kfold: Number of folds for cross-validation.
        gpu_id: GPU ID to use for training (0 for first GPU).
        label_weight: Weight for each label in loss calculation.
        
    Raises:
        RuntimeError: If training fails.
    """
    if label_weight is None:
        label_weight = [0.3, 0.7]
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Starting training with {kfold}-fold cross-validation...")
    print(f"Training parameters:")
    print(f"  Epochs: {epochs}")
    print(f"  Batch size: {batch_size}")
    print(f"  Learning rate: {learning_rate}")
    print(f"  Early stopping: {early_stopping}")
    print(f"  Label weights: {label_weight}")
    print(f"  Output directory: {output_dir}")
    
    try:
        clf = MolTrain(
            task="multilabel_regression",
            data_type="molecule_all_h",
            epochs=epochs,
            batch_size=batch_size,
            metrics=["mse"],
            learning_rate=learning_rate,
            early_stopping=early_stopping,
            save_path=str(output_dir),
            kfold=kfold,
            gpu_id=gpu_id,
            loss_key="weighted_mse",
            label_weight=label_weight,
        )
        
        clf.fit(data=data)
        print(f"Training completed successfully. Model saved to: {output_dir}")
        
    except Exception as e:
        raise RuntimeError(f"Training failed: {e}") from e


def main():
    """Main function to run model training."""
    parser = argparse.ArgumentParser(
        description="Train UniMol model on peptide data"
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="Path to training data pickle file (default: ../data/23pep.pkl)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output directory for trained model (default: ../output)"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=500,
        help="Number of training epochs (default: 500)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
        help="Batch size for training (default: 128)"
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=1e-4,
        help="Learning rate (default: 1e-4)"
    )
    parser.add_argument(
        "--early-stopping",
        type=int,
        default=20,
        help="Early stopping patience (default: 20)"
    )
    parser.add_argument(
        "--kfold",
        type=int,
        default=4,
        help="Number of folds for cross-validation (default: 4)"
    )
    parser.add_argument(
        "--gpu-id",
        type=int,
        default=0,
        help="GPU ID to use (default: 0)"
    )
    parser.add_argument(
        "--label-weight",
        type=float,
        nargs=2,
        default=[0.3, 0.7],
        help="Weights for labels in loss calculation (default: 0.3 0.7)"
    )
    
    args = parser.parse_args()
    
    # Setup paths
    project_root = get_project_root()
    data_dir = project_root / "data"
    
    data_path = Path(args.data) if args.data else data_dir / "23pep.pkl"
    output_dir = Path(args.output) if args.output else project_root / "output"
    
    try:
        # Load training data
        print(f"Loading training data from: {data_path}")
        data = load_training_data(data_path)
        
        # Train model
        train_model(
            data=data,
            output_dir=output_dir,
            epochs=args.epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate,
            early_stopping=args.early_stopping,
            kfold=args.kfold,
            gpu_id=args.gpu_id,
            label_weight=args.label_weight
        )
        
    except (FileNotFoundError, RuntimeError, pickle.UnpicklingError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()