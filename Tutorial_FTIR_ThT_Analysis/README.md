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
  Usage: `python fit_impr_plot.py example/CITI.csv`
* `clustering_spectral_lin_all.ipynb` — Spectral clustering (RBF kernel) of
  peptides in 2D (f_beta from FTIR, log10(F/F0) from the ThT assay) space to
  assign binary amyloid-forming / non-amyloid ground-truth labels, and
  reports the clustering silhouette score. Reads `sum_complie_all.csv` and
  `tht_60.csv` (both included in this folder).
* `sum_complie_all.csv` — compiled per-peptide results for the 60
  characterized peptides: `beta` = f_beta (%) from FTIR deconvolution,
  `para` = parallel beta-sheet content (%) derived from the FTIR high/low
  frequency component ratio, `Amy` = binary amyloid-forming ground-truth
  label (1/0) from spectral clustering of the (beta, ThT) data.
* `tht_60.csv` — Thioflavin-T (ThT) fluorescence intensities (`f1`-`f3` =
  triplicate raw measurements, `sem` = standard error of the mean) and the
  derived `logF/F0` (relative fluorescence enhancement vs. the ThT-only
  blank) with its propagated SEM (`sem_log`).
* `example/` — one complete worked example (peptide CITI): the raw FTIR
  input spectrum (`CITI.csv`), the fitted component curves
  (`CITI_fit_FTIR_SourceData.csv`), and the resulting figure
  (`CITI_fit.png`), i.e. the input/output pair for `fit_impr_plot.py`.

## Prerequisites

```bash
pip install numpy pandas matplotlib scipy scikit-learn adjustText
```

## Usage

1. **FTIR deconvolution** (per peptide): run the script on a raw FTIR
   spectrum (`<SEQ>.csv`, columns `wavenumber,Absorbance`), e.g. the bundled
   example:

```bash
python fit_impr_plot.py example/CITI.csv
```

   This produces `<SEQ>_fit.png` and `<SEQ>_fit_FTIR_SourceData.csv`
   (fitted component curves) alongside the input file, matching the files
   already provided in `example/` for CITI.

2. **Spectral clustering**: run all cells of
   `clustering_spectral_lin_all.ipynb` (it reads `sum_complie_all.csv` and
   `tht_60.csv` from this same folder). This reproduces the
   amyloid-forming/non-amyloid classification (Figs. 3g and 4c of the
   manuscript) and reports the silhouette score.

## Related resources

* PACE simulation results (AP/ALT scores, final assembly structures) and
  synthesis Certificates of Analysis: Zenodo,
  https://doi.org/10.5281/zenodo.21253232
* PACE MD simulation automation: [`Tutorial_Automated_Peptide_MD_PACE`](../Tutorial_Automated_Peptide_MD_PACE)
* Uni-Mol machine learning model: [`Self_Assembling_Amyloid_Like_Peptides_Predictor`](../Self_Assembling_Amyloid_Like_Peptides_Predictor)
* Raw full FTIR/ThT spectra (beyond the CITI example above) and SEM/XRD
  imaging data for all characterized peptides are available from the
  corresponding authors upon reasonable request.
