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

> **TL;DR**: Don't forget Git LFS (Step 0a) or you'll silently get fake/empty model files. If you have a Linux machine with an NVIDIA GPU, CUDA 11.7 or 11.8, and are willing to use Python 3.8/3.9/3.10 + PyTorch 2.0.0, you can install everything with pre-compiled wheels in ~6 commands (Step 0-5 below, Option A for Uni-Core). Otherwise (different CUDA/PyTorch version, or no GPU), you'll need to build Uni-Core from source (Option B) — it's more steps but still straightforward if you follow them in order.

### Prerequisites

- **OS**: Linux only. Uni-Core ships custom CUDA kernels that are only pre-built/tested for Linux; macOS/Windows are not supported (WSL2 on Windows should work like Linux).
- **Python**: 3.7 – 3.10 (3.8/3.9/3.10 recommended, since those are the versions with pre-compiled Uni-Core wheels available).
- **GPU**: An NVIDIA GPU + driver is strongly recommended. CPU-only is possible (it's actually the default when building Uni-Core from source, see Step 3) but will be much slower.
- **Disk/network**: You'll clone/download a few packages (Uni-Core, this repo), so make sure you have normal internet access (or a local mirror) for `pip`/`git`.

### Step 0a: Make sure Git LFS is set up before you clone (easy to miss!)

This repository stores large binary files — the pre-trained weights under `trained_model/`/`validated_model/` (`*.pth`, ~190 MB each) and `data/seq2smi_1to5.json` (~400 MB) — via **Git LFS**, not as plain Git objects. If `git-lfs` isn't installed when you clone, Git silently checks out tiny (~130-byte) *pointer files* instead of the real content, e.g.:

```
version https://git-lfs.github.com/spec/v1
oid sha256:8ab7211efb5b2bccb4429b3a228caf6bfa32a579edb058cf7ff0bbd1652af71f
size 190447415
```

No error is shown at clone time — the failure only shows up later as a confusing `Error: Expecting value: line 1 column 1 (char 0)` (a JSON-parse error) when `predict.py` tries to read one of these pointer files as if it were real data/weights. To avoid this, install Git LFS **before** cloning (or run the pull step below right after cloning if you already have a broken checkout):

```bash
# macOS
brew install git-lfs
# Debian/Ubuntu
sudo apt-get install git-lfs

git lfs install                 # one-time, per machine
git clone https://github.com/hanlab-computChem/hanlab.git
cd hanlab
git lfs pull --include "Self_Assembling_Amyloid_Like_Peptides_Predictor/*"   # fetch the real files for this subproject
```

If you already have a checkout with placeholder files (check with `ls -la trained_model/*.pth` — real files are ~190 MB, pointer files are ~130 bytes), just run `git lfs install && git lfs pull --include "Self_Assembling_Amyloid_Like_Peptides_Predictor/*"` from the repo root to fix it in place. Note this downloads roughly 1-2 GB and can take a while depending on your connection.

### Step 0: Create an isolated environment (recommended)

Using a dedicated conda environment avoids version conflicts with other projects:

```bash
conda create -n unimol_peptide python=3.9 -y
conda activate unimol_peptide
```

(Any Python 3.7–3.10 works; 3.9 is used here because it has a matching pre-compiled Uni-Core wheel — see Step 3.)

### Step 1: Check your CUDA version

The CUDA version used to build PyTorch and the CUDA version used to build Uni-Core **must match**, or you'll get import/runtime errors later. Check what's available on your machine first:

```bash
nvidia-smi        # top-right corner shows the max CUDA version your driver supports
nvcc --version     # shows the CUDA toolkit version actually installed (may differ from nvidia-smi)
```

- If you don't have `nvcc`/CUDA toolkit installed system-wide, that's fine — PyTorch ships its own CUDA runtime, so you only need `nvidia-smi` to show a version >= the CUDA version you plan to install PyTorch with.
- Pre-compiled Uni-Core wheels are currently only available for **CUDA 11.7 / 11.8 + PyTorch 2.0.0**. If your driver supports a newer CUDA (e.g. 12.x), it can usually still run CUDA 11.8 binaries (CUDA is backward compatible), so prefer installing PyTorch 2.0.0 + cu118 to stay on the easy path (Option A in Step 3).

### Step 2: Install PyTorch

Install PyTorch **before** the rest of the dependencies, matching the CUDA version you chose above. For example, for CUDA 11.8:

```bash
pip install torch==2.0.0 --index-url https://download.pytorch.org/whl/cu118
```

For CUDA 11.7:

```bash
pip install torch==2.0.0 --index-url https://download.pytorch.org/whl/cu117
```

If you need a different PyTorch/CUDA combination (e.g. no pre-compiled Uni-Core wheel matches your setup), install PyTorch per the [official PyTorch install guide](https://pytorch.org/get-started/locally/) instead, then build Uni-Core from source in Step 3 (Option B).

Verify the install and that CUDA is visible to PyTorch:

```bash
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```

You should see something like `2.0.0 11.8 True`. If `torch.cuda.is_available()` prints `False` on a GPU machine, stop here and fix your driver/CUDA setup before continuing — nothing downstream will work correctly otherwise.

### Step 3: Install Uni-Core

This is the step most people get stuck on, so read carefully. [Uni-Core](https://github.com/dptech-corp/Uni-Core) provides fused CUDA kernels that Uni-Mol (and therefore this package) depends on.

**Option A — pre-compiled wheel (fastest, recommended if it matches your setup)**

Pre-compiled wheels are only published for **Python 3.8/3.9/3.10 + PyTorch 2.0.0 + CUDA 11.7/11.8**. Pick the file matching your Python/CUDA combo from the [Uni-Core Releases page](https://github.com/dptech-corp/Uni-Core/releases) and `pip install` it directly by URL, e.g. for Python 3.9 + CUDA 11.8:

```bash
pip install https://github.com/dptech-corp/Uni-Core/releases/download/0.0.3/unicore-0.0.1+cu118torch2.0.0-cp39-cp39-linux_x86_64.whl
```

The filename encodes the requirements: `unicore-<unicore_version>+cu<CUDA_version>torch<torch_version>-cp<python_version>-...`. Match every field (`cp39` = Python 3.9, `cu118` = CUDA 11.8, `torch2.0.0` = PyTorch 2.0.0) against your environment from Steps 1-2, or the import will fail with an ABI/symbol mismatch.

If none of the published wheels match your Python/PyTorch/CUDA combination, use Option B instead.

**Option B — build from source (use this if your CUDA/PyTorch/Python version isn't covered by a wheel, or you're on CPU-only)**

Building requires a C++ compiler and (for the CUDA extension) the CUDA toolkit matching your PyTorch build:

```bash
# Make sure these are on PATH/set correctly before building:
#   gcc/g++ (a version compatible with your CUDA toolkit, e.g. gcc 9-11 for CUDA 11.x)
#   CUDA_HOME pointing at your CUDA toolkit install, e.g.:
export CUDA_HOME=/usr/local/cuda-11.8   # adjust to match your installed toolkit
export PATH=$CUDA_HOME/bin:$PATH

git clone https://github.com/dptech-corp/Uni-Core.git
cd Uni-Core
python setup.py install --enable-cuda-ext
cd ..
```

> **Verified 2026-07 against the current Uni-Core `main` branch**: the CUDA extension flag has flipped since the 0.0.3 release (which is what the pre-compiled wheels above were built from). On `main`, the CUDA extension is now **disabled by default**, and you must pass `--enable-cuda-ext` explicitly to build it — the old `--disable-cuda-ext` flag from the upstream Uni-Core README no longer exists and will fail with `error: option --disable-cuda-ext not recognized`. If you hit that error, it means the flag has changed again; try building with no flag at all (CPU-only) or check `Uni-Core/setup.py` directly for the current flag name.

On a machine with no GPU/CUDA toolkit at all (CPU-only), just omit the flag — this is now the default:

```bash
python setup.py install
```

This is expected to print several `<kernel_name> is not installed corrected` lines on first `import unicore` afterwards (e.g. `fused_layer_norm is not installed corrected`) — that's normal for a CPU-only/no-CUDA-ext build, not a failure; see Troubleshooting.

Building from source typically takes a few minutes for a CPU-only build (longer with `--enable-cuda-ext`, since it compiles fused CUDA kernels). If it fails, see Troubleshooting below.

**Verify Uni-Core installed correctly:**

```bash
python -c "import unicore; print('unicore OK')"
```

### Step 4: Install this package's remaining Python dependencies

```bash
cd Self_Assembling_Amyloid_Like_Peptides_Predictor
pip install -r requirements.txt
```

This installs `rdkit`, `pandas`, `transformers`, and the other libraries listed in `requirements.txt` (Uni-Core is intentionally *not* in this file since it must be installed separately per Step 3).

### Step 5: Install this package itself

```bash
python setup.py install
# or, for an editable install that picks up local code changes without reinstalling:
# pip install -e .
```

### Step 6 (optional): Install `confgen` for higher-quality 3D conformers

The prediction functionality can optionally use the `confgen` command-line tool — the CONFORGE conformer generator bundled with [CDPKit](https://cdpkit.org/) ([molinfo-vienna/CDPKit](https://github.com/molinfo-vienna/CDPKit)) — to generate 3D conformations.

- `confgen` is a CLI binary, not just the `pip install cdpkit` Python bindings, so download it from the [CDPKit installer packages](https://github.com/molinfo-vienna/CDPKit/releases) (macOS/Linux/Windows) or build CDPKit from source, then add its `Bin` folder to your `PATH` (see the [CDPKit installation docs](https://cdpkit.org/installation.html)).
- Check it's on your `PATH`: `confgen --help` should print usage info.
- **This step is optional.** If `confgen` isn't found in `PATH`, `examples/predict.py` automatically falls back to an RDKit-based pipeline (`Chem.AddHs` + `EmbedMultipleConfs` + MMFF optimization), so everything still works as long as RDKit is available. Either way, conformers are generated with explicit hydrogens first for accurate 3D geometry, and the hydrogens are then stripped before the coordinates reach the model (see `remove_hs` under Technical Details).

### Step 7: Verify the full installation

Run a quick end-to-end smoke test using the bundled pre-trained model:

```bash
cd examples
python predict.py -p FFFF WWWW YYYY -o predict_result.csv
cat predict_result.csv
```

If installation succeeded, you should see output similar to:

```csv
Peptide,Pred_AP,Pred_SHB
FFFF,3.8230,0.5959
WWWW,3.1570,0.2713
YYYY,3.2162,0.2086
```

If this runs without errors and produces a CSV with numeric predictions, the installation is complete.

### Troubleshooting

The entries below marked ✅ were reproduced and confirmed while writing this guide (in a clean Python 3.9 + PyTorch 2.0.0 CPU-only container); the rest are common failure modes documented from the underlying tools' own issue trackers.

- ✅ **`Error: Expecting value: line 1 column 1 (char 0)` while running `predict.py`/`train.py`**: You have Git LFS *pointer files* instead of real data/weights (see Step 0a) — likely because `git-lfs` wasn't installed when you cloned. Run `git lfs install && git lfs pull --include "Self_Assembling_Amyloid_Like_Peptides_Predictor/*"` from the repo root and re-run. Quick sanity check: `ls -la trained_model/*.pth` should show ~190 MB files, not ~130 bytes.
- ✅ **`ModuleNotFoundError: No module named 'sklearn'`**: `scikit-learn` (imported as `sklearn`) and `joblib` are required by `unimol_tools` but were missing from `requirements.txt` in earlier versions of this doc/file — they've since been added. If you still hit this (e.g. you installed dependencies manually rather than via `requirements.txt`), just run `pip install scikit-learn joblib`.
- ✅ **`NameError: name 'LRScheduler' is not defined`** (raised deep inside `transformers/trainer_pt_utils.py` when importing `unimol_tools`): You installed an unpinned/too-new `transformers` that assumes `torch>=2.1`'s `torch.optim.lr_scheduler.LRScheduler`, which doesn't exist in `torch==2.0.0`. Fix with `pip install "transformers==4.30.2"` (already pinned in `requirements.txt`).
- ✅ **`<kernel_name> is not installed corrected` printed when you `import unicore`** (e.g. `fused_layer_norm is not installed corrected`): This is expected and harmless if you built/installed Uni-Core without the CUDA extension (CPU-only, or you didn't pass `--enable-cuda-ext`) — those fused kernels just aren't available, and pure-PyTorch fallbacks are used instead. It is not something you need to "fix" unless you specifically need the fused-kernel speedup on GPU.
- **`error: option --disable-cuda-ext not recognized` when building Uni-Core from source**: You're on a newer Uni-Core `main` checkout where the flag was renamed/flipped — use `--enable-cuda-ext` instead (or nothing, for CPU-only). See the note in Step 3, Option B.
- **`ImportError` / `undefined symbol` when importing `unicore`**: Your Uni-Core wheel/build doesn't match your installed PyTorch/CUDA/Python version. Re-check `torch.__version__`, `torch.version.cuda`, and `python --version`, then re-install Uni-Core with the exact matching combination (Step 3).
- **`torch.cuda.is_available()` returns `False`**: Usually a driver/CUDA mismatch, not a Uni-Core problem. Fix this before installing Uni-Core — verify with `nvidia-smi` and reinstall the correct PyTorch CUDA build (Step 2).
- **Building Uni-Core from source fails with a CUDA/nvcc error**: Make sure `CUDA_HOME` points at a CUDA toolkit whose version matches `torch.version.cuda`, and that `nvcc --version` under `$CUDA_HOME/bin` reports the same major.minor version. Mismatched toolkit versions are the most common cause of build failures.
- **Building Uni-Core from source fails with a compiler/`gcc` error**: Your `gcc`/`g++` version may be too new or too old for your CUDA toolkit (e.g. CUDA 11.x generally needs gcc <= 11). Install a compatible compiler version (e.g. via `conda install -c conda-forge gcc=9 gxx=9`) and retry.
- **No GPU available**: Build Uni-Core with plain `python setup.py install` (no `--enable-cuda-ext`, Step 3 Option B) and install a CPU-only PyTorch build. Predictions/training will still run, just on CPU.
- **`confgen: command not found`**: This is expected if you skipped Step 6 — it's optional and the code falls back to RDKit automatically. Only fix this if you specifically want CDPKit-quality conformers.
- **`pip install -r requirements.txt` seems to hang or is very slow**: This is usually a slow PyPI mirror for large packages like `transformers`/`wandb`. Consider using a local/regional PyPI mirror.

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

Note that `examples/train.py` doesn't currently expose a `--remove-hs` flag and relies on `MolTrain`'s default, so if you're retraining to reproduce the released models, call `MolTrain` through the Python API instead and pass `remove_hs=True` explicitly (see below).

## Using Python API

### Prediction

```python
from unimol_tools import MolPredict
from pathlib import Path

# Load pre-trained model
model_path = Path("trained_model")
clf = MolPredict(load_model=str(model_path))

# Prepare data (dictionary format). Generate atoms/coordinates with
# hydrogens included (e.g. via RDKit's AddHs + conformer embedding) for
# accurate geometry — MolPredict will strip them automatically per the
# loaded model's remove_hs setting.
data = {
    'atoms': [['C', 'C', 'O', 'N', 'H', 'H', ...], ...],  # List of atom symbols
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
    remove_hs=True,  # strip hydrogens after conformer generation, to match the released models
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

- **Base Model**: Uni-Mol (Universal 3D Molecular Representation Learning Framework), fine-tuned from the all-hydrogen pretrained checkpoint (`data_type=molecule_all_h`)
- **Task Type**: Multilabel regression (multilabel_regression)
- **Data Format**: Molecular 3D coordinates generated with explicit hydrogens (for accurate conformer geometry), then with hydrogens removed before tokenization — i.e. `remove_hs=True`, as set in each model's `config.yaml`. `MolPredict` picks this up automatically from the loaded model directory, so no extra step is needed for prediction; when training your own model with `MolTrain`, pass `remove_hs=True` explicitly to match it
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
