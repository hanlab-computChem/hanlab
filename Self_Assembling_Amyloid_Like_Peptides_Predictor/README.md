# Self-Assembling Amyloid-Like Peptides Predictor

A property prediction tool for self-assembling amyloid-like peptides based on the Uni-Mol framework. This tool uses deep learning models to predict two key properties of short peptide sequences:
- **AP (Amyloid Propensity)**: Amyloid-forming propensity
- **SHB (Self-Assembly/β-sheet propensity)**: Self-assembly and β-sheet formation tendency

## Features

- **Multilabel Regression Prediction**: Simultaneously predict both AP and SHB properties of peptide sequences
- **Pre-trained Models**: Pre-trained Uni-Mol-based models ready for direct use
- **Model Training**: Support for training new models with custom datasets
- **Cross-Validation**: K-fold cross-validation support during training (default: 4-fold)
- **3D Conformation Generation**: Automatically convert peptide sequences to SMILES and generate 3D coordinates

## Installation

### Prerequisites

- OS: Linux (Uni-Core ships custom CUDA kernels that are only pre-built for Linux; other OSes are not officially supported)
- Python 3.7 – 3.10
- PyTorch >= 2.0.0, with a CUDA build that matches the CUDA toolkit used to build Uni-Core (see below)
- RDKit (used for SMILES → 3D conformer generation, `rdkit>=2022.9.3`)
- [Uni-Core](https://github.com/dptech-corp/Uni-Core) - Please install Uni-Core first

### Installation Steps

1. **Install PyTorch** matching your CUDA version first (see the [PyTorch install guide](https://pytorch.org/get-started/locally/)), then **install the remaining dependencies**:
```bash
pip install -r requirements.txt
```

2. **Install Uni-Core** (if not already installed). Pick one of the following, and make sure the CUDA version matches the one used to build your PyTorch:
```bash
# Option A: pre-compiled wheel (fastest, see the Releases page for the
# wheel matching your python/torch/cuda version)
# https://github.com/dptech-corp/Uni-Core/releases
pip install unicore-<version>+cu<cuda>torch<torch_version>-<abi_tag>-linux_x86_64.whl

# Option B: build from source
git clone https://github.com/dptech-corp/Uni-Core.git
cd Uni-Core
python setup.py install
# use `python setup.py install --disable-cuda-ext` on machines without a GPU/CUDA toolchain
```

3. **Install this package**:
```bash
python setup.py install
```

### Additional Dependencies

The prediction functionality can optionally use the `confgen` tool (DP Technology's internal conformer generator) to generate 3D conformations. If `confgen` is not found in your system `PATH`, `examples/predict.py` automatically falls back to an RDKit-based conformer generation pipeline (`Chem.AddHs` + `EmbedMultipleConfs` + MMFF optimization), so installing `confgen` is optional as long as RDKit is available.

### Hydrogen Handling (`remove_hs`) — Required for Reproducible Predictions

Both released models (`trained_model/` and `validated_model/`) were trained with `data_type="molecule_all_h"` and `remove_hs=False` (see `config.yaml` in each model directory). This means **all hydrogen atoms must be kept, not removed**, in the 3D structures fed to the model:

- The bundled `examples/predict.py` and `examples/train.py` scripts already do the right thing by default: `Chem.AddHs`/`confgen` always produce explicit hydrogens, and `remove_hs` is never set to `True`, so no action is needed to reproduce the documented results when using these scripts as-is.
- If you build your own `atoms`/`coordinates` input for the Python API (see below), make sure hydrogen atoms are included and **not** stripped out — omitting or removing hydrogens will silently change the input representation and will not reproduce the reported AP/SHB predictions.
- If you retrain or fine-tune a model, keep `remove_hs=False` (the default in `MolTrain`/`MolPredict`) unless you intentionally want a no-hydrogen (`molecule_no_h`) model, in which case predictions must also be made with `remove_hs=True` to stay consistent.

## Quick Start

### Prediction with Pre-trained Models

Use the `examples/predict.py` script to predict peptide sequences:

```bash
cd examples
python predict.py -p FFFF WWWW YYYY -o predict_result.csv
```

**Parameters**:
- `-p, --peptides`: List of peptide sequences to predict (space-separated)
- `-o, --output`: Output CSV file path (optional, defaults to `predict_result.csv`)

**Output Format**:
```csv
Peptide,Pred_AP,Pred_SHB
FFFF,3.8230,0.5959
WWWW,3.1570,0.2713
YYYY,3.2162,0.2086
```

**Notes**:
- Currently supports peptide sequences of 1-5 amino acids
- Peptide sequences are automatically converted to uppercase
- Sequences not in the supported mapping table will be skipped with a warning

### Training New Models

Use the `examples/train.py` script to train new models:

```bash
cd examples
python train.py --data ../data/3pep_assembly_train.pkl --output ../output --epochs 500 --batch-size 128
```

**Main Parameters**:
- `--data`: Path to training data pickle file (default: `../data/23pep.pkl`)
- `--output`: Model output directory (default: `../output`)
- `--epochs`: Number of training epochs (default: 500)
- `--batch-size`: Batch size (default: 128)
- `--learning-rate`: Learning rate (default: 1e-4)
- `--kfold`: Number of cross-validation folds (default: 4)
- `--gpu-id`: GPU ID to use (default: 0)
- `--label-weight`: Label weights for weighted loss function (default: 0.3 0.7)

**Training Data Format**:
Training data should be a dictionary or pickle file containing the following fields:
- `atoms`: List of lists of atom symbols
- `coordinates`: List of 3D coordinate arrays
- `target`: Label array containing both AP and SHB values

## Using Python API

### Prediction

```python
from unimol_tools import MolPredict
from pathlib import Path

# Load pre-trained model
model_path = Path("trained_model")
clf = MolPredict(load_model=str(model_path))

# Prepare data (dictionary format)
# NOTE: 'atoms'/'coordinates' must include explicit hydrogen atoms
# (do NOT remove hydrogens) to match how the released models were
# trained (remove_hs=False, data_type='molecule_all_h'). See
# "Hydrogen Handling" above.
data = {
    'atoms': [['C', 'C', 'O', 'N', 'H', 'H', ...], ...],  # List of atom symbols, including H
    'coordinates': [np.array([...]), ...],      # 3D coordinate arrays
}

# Make predictions
predictions = clf.predict(data=data)
# predictions shape: (n_samples, 2) - [AP, SHB]
```

### Training

```python
from unimol_tools import MolTrain

# Initialize trainer
clf = MolTrain(
    task='multilabel_regression',
    data_type='molecule_all_h',
    remove_hs=False,  # keep all hydrogens, required to reproduce the released models
    epochs=500,
    batch_size=128,
    metrics=['mse'],
    learning_rate=1e-4,
    early_stopping=20,
    kfold=4,
    gpu_id=0,
    loss_key='weighted_mse',
    label_weight=[0.3, 0.7],  # Weights for AP and SHB
    save_path='./output'
)

# Train model
clf.fit(data=training_data)
```

## Project Structure

```
Self_Assembling_Amyloid_Like_Peptides_Predictor/
├── examples/              # Example scripts
│   ├── train.py          # Training script
│   └── predict.py        # Prediction script
├── unimol_tools/         # Core toolkit
│   ├── config/           # Configuration files
│   ├── data/             # Data processing modules
│   ├── models/           # Model definitions
│   ├── tasks/            # Training tasks
│   ├── utils/            # Utility functions
│   ├── train.py          # Training interface
│   ├── predict.py        # Prediction interface
│   └── predictor.py      # Predictor implementation
├── data/                 # Data files
│   └── seq2smi_1to5.json # Peptide sequence to SMILES mapping
├── trained_model/        # Fine-tuned models used for production prediction such as tetrapeptide prediction
│   ├── config.yaml       # Model configuration
│   ├── model_*.pth       # Model weight files
│   └── target_scaler.ss  # Target value scaler
├── validated_model/      # Fine-tuned models used for model validation via tripeptide prediction
│   ├── config.yaml       # Model configuration
│   ├── model_*.pth       # Model weight files
│   └── target_scaler.ss  # Target value scaler
├── requirements.txt      # Python dependencies
├── setup.py             # Installation script
└── README.md            # This document
```

## Technical Details

### Model Architecture

- **Base Model**: Uni-Mol (Universal 3D Molecular Representation Learning Framework)
- **Task Type**: Multilabel regression (multilabel_regression)
- **Data Format**: Molecular 3D coordinates with hydrogens (molecule_all_h)
- **Hydrogen Handling**: `remove_hs=False` — hydrogens are kept, not removed (see "Hydrogen Handling" note in the Installation section above)
- **Loss Function**: Weighted mean squared error (weighted_mse)
- **Evaluation Metric**: Mean squared error (MSE)

### Training Configuration

- **Cross-Validation**: 4-fold random split
- **Batch Size**: 128
- **Learning Rate**: 1e-4
- **Max Epochs**: 500
- **Early Stopping Patience**: 20 epochs
- **Label Weights**: AP=0.3, SHB=0.7

## Related Resources

- **Uni-Mol Documentation**: https://unimol.readthedocs.io/en/latest/
- **Uni-Mol Paper**: [Uni-Mol: A Universal 3D Molecular Representation Learning Framework](https://openreview.net/forum?id=6K2RM6wVqKu)
- **Uni-Core Repository**: https://github.com/dptech-corp/Uni-Core
