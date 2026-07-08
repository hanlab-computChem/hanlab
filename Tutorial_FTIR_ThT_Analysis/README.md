# Tutorial: FTIR Deconvolution and ThT/FTIR Spectral Clustering for Amyloid-Forming Peptides

This folder contains the scripts used to (i) deconvolute ATR-FTIR spectra of
short peptides into secondary-structure components and (ii) classify
peptides into amyloid-forming/non-amyloid groups from their combined ThT and
FTIR results, as described in the Methods section ("Attenuated total
reflectance-Fourier transform infrared spectroscopy" and "Clustering
analysis") of the associated manuscript.

## Contents

* `fit_impr_plot.py` — Gaussian deconvolution of the Amide I region
  (1580-1720 cm-1) of an ATR-FTIR spectrum into 7 component peaks (baseline,
  cross-beta, random coil, helix/turn, anti-parallel sheet,
  primary-amide/side-chain, and pre-proline/C=O stretch contributions).
  Computes the beta-sheet fraction f_beta and the parallel/anti-parallel
  ratio R_a/p, and generates a fit plot.
* `clustering_spectral_lin_all.ipynb` — Spectral clustering (RBF kernel) of
  peptides in 2D (f_beta from FTIR, log10(F/F0) from the ThT assay) space to
  assign binary amyloid-forming / non-amyloid ground-truth labels, and
  reports the clustering silhouette score.

## Input data

These scripts consume the compiled experimental data (raw ATR-FTIR spectra,
compiled FTIR/ThT summary tables, and a worked example) archived at Zenodo:
**https://doi.org/10.5281/zenodo.21253232** (folder `04_experimental_data/`
in that record).

## Prerequisites

```bash
pip install numpy pandas matplotlib scipy scikit-learn adjustText
```

## Usage

1. **FTIR deconvolution** (per peptide): download the raw FTIR spectrum
   (`<SEQ>.csv`, columns `wavenumber,Absorbance`) from the Zenodo
   `04_experimental_data/example/` folder (or your own ATR-FTIR export in
   the same format), then run:

```bash
python fit_impr_plot.py example/<SEQ>.csv
```

   This produces `<SEQ>_fit.png` and `<SEQ>_fit_FTIR_SourceData.csv`
   (fitted component curves) alongside the input file.

2. **Spectral clustering**: place `sum_complie_all.csv` and `tht_60.csv`
   (from the Zenodo `04_experimental_data/` folder) in the same directory as
   the notebook, then run all cells of `clustering_spectral_lin_all.ipynb`.
   This reproduces the amyloid-forming/non-amyloid classification (Figs. 3g
   and 4c of the manuscript) and reports the silhouette score.

## Related resources

* Compiled experimental data (input files for the scripts above), synthesis
  Certificates of Analysis, and PACE simulation results: Zenodo,
  https://doi.org/10.5281/zenodo.21253232
* PACE MD simulation automation: [`Tutorial_Automated_Peptide_MD_PACE`](../Tutorial_Automated_Peptide_MD_PACE)
* Uni-Mol machine learning model: [`Self_Assembling_Amyloid_Like_Peptides_Predictor`](../Self_Assembling_Amyloid_Like_Peptides_Predictor)
