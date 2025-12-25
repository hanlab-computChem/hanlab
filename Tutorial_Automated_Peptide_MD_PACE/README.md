# Tutorial: Automated GROMACS Simulation Setup for Large Sequence Batches

This guide outlines the process for automatically generating GROMACS simulation files for multiple peptide sequences using the **PACE Force Field** and the **Auto-builder** tool.

## Prerequisites

Before starting, ensure your environment meets the following requirements:

* **GROMACS:** Version 5.0 or later.
* **Compiler:** C compiler (GCC/GNU preferred).
* **Python:** Version 3.x.
* **Python Libraries:** `PeptideConstructor`, `biopython`, and `PeptideBuilder`.
* **Installation:** `pip install PeptideConstructor biopython PeptideBuilder`

## I. PACE Force Field Installation

1. **Download:** Get the PACE package from the [Han Lab GitHub repository](https://github.com/hanlab-computChem/hanlab/blob/master/Tutorial_PACE-ASM/Tutorial_PACE-ASM.zip).


2. **Extract:** Unzip `Tutorial_PACE-ASM.zip` and navigate into the `Tutorial_PACE-ASM/` directory.


3. **Setup:** Run the setup script:
```bash
bash setup.sh
```


This generates a script named `prepare.sh`, which is required to construct GROMACS topology files for PACE.



> 
> **Note:** A detailed guide is included in this package for advanced troubleshooting.
> 
> 

---

## II. Auto-builder Installation

1. **Download:** Download the [automated building tool](https://github.com/hanlab-computChem/hanlab/blob/master/Tutorial_Automated_Peptide_MD_PACE/Auto_build).


2. **Extract:** Unzip the package and enter the `Auto_build/` directory. You will see two subdirectories: `build/` and `setup/`.


3. **Integration:** Copy the `prepare.sh` script generated in Section I into the `build/` directory (assuming PACE and Auto-building tool are in the same parent directory):


```bash
cp ../Tutorial_PACE-ASM/prepare.sh ./build/
```






---

## III. Automatic Generation of Simulation Files

1. **Environment Check:** Ensure your OS environment variables are correctly set so that `gmx` (GROMACS) and `python3` are accessible from the command line.


2. **Prepare Sequences:** Navigate to the `build/` directory and create a text file (e.g., `seq.txt`) containing the peptide sequences to be simulated (one per line).


**Example `seq.txt`:**
```text
AAA
AASF
```





3. **Run the Script:** Execute the batch builder in the background:
```bash
bash build_byFile.sh seq.txt > log 2>&1 &
```






---

## IV. Troubleshooting

### 1. Verification Test

To verify the installation, run a test with a single sequence (e.g., `AAA`). Check the generated `log` file.

A successful run should include:

* **Structure Generation:** `==> Peptide 1A-2A-3A has been generated and saved in AAA.pdb` 


* **Topology Generation:**
```text
--- PLEASE NOTE ---
You have successfully generated a topology from: AAA-cap.pdb.
The Pace-asm force field and the cgWater water model are used.
--- ETON ESAELP ---
```






### 2. Common Issues

If the topology fails or simulations do not start:

* **Topology building Errors:** Consult the PACE documentation


* **Environment Errors:** Double-check that your GROMACS and Python paths are correctly exported.


* **Version Mismatch:** Minor syntax differences between GROMACS versions can cause errors. If the simulation fails during pre-equilibrium, inspect and edit the configuration templates in: `setup/gentop_homo_num.sh`.



---

## V. Results

Upon successful completion, all simulation files (`.tpr`) will be moved to the following directory for production runs: `setup/runtpr/`.

---

## VI. Automated Analysis

You can automate the calculation of structural properties for your assemblies using the `Auto_analysis` tool. This tool processes batches of simulation results simultaneously.

### 1. Preparation

1. **Download:** Get the tool via [link](https://github.com/hanlab-computChem/hanlab/blob/master/Tutorial_Automated_Peptide_MD_PACE/Auto_analysis).
2. **Extract:** Unzip the file and enter the `Auto_analysis/` directory.
3. **Environment:** Ensure your environment variables (GROMACS and Python) are identical to those used in the simulation setup.
4. **Dependencies:** Install the MDAnalysis library:
```bash
pip install mdanalysis

```






### 2. Data Organization

For the script to function correctly, your simulation data must be organized in a specific directory structure. Each sequence should have its own folder containing the `.tpr` and `.xtc` files:

* `PATH_TO_SIM_DATA/AAA/` (contains `AAA.tpr`, `AAA.xtc`) 


* `PATH_TO_SIM_DATA/AAD/` (contains `AAD.tpr`, `AAD.xtc`) 


* `PATH_TO_SIM_DATA/AAE/` (contains `AAE.tpr`, `AAE.xtc`) 



### 3. Running the Analysis

1. **Configure Parameters:** Open `ana_all.sh` and locate the parameter section as follows.

```txt
##MAJOR parameters for analysis##
# length of peptide sequence, e.g., 4 for tetrapeptides seqLen=3
# frame frequency for analysis jump_step=10
# select correct directories where the results for each sequence
# directory where the simulation results are stored
# This is for demonstration
RTDIR=$WORKDIR/example
# for acutal application, please use the following statement
#RTDIR=PATH_TO_SIM_RESULT
# directory where the analysis results are stored
# This is for demonstration
TGTDIR=$WORKDIR/example/Result
# for acutal application, please use the following statement
#TGTDIR=PATH_TO_ANA_RESULT
# name of GMX executable usually gmx or gmx_mpi
gmx_mpi=gmx
END OF parameter SESSION
```

Set your analysis parameters according to the comments provided within the script.



2. **Execute:** Run the analysis script:

```bash
bash ana_all.sh

```



3. **Verification:** To confirm everything processed correctly, run the script a second time:

```bash
bash ana_all.sh
```

If successful, the output should contain no "missing" warnings.

### 4. Results

The analysis outputs (plots, data files, etc.) will be generated in a corresponding directory structure: `PATH_TO_RESULTS/AAA/`, `PATH_TO_RESULTS/AAD/`, etc..

> 
> **Note:** A reference example of the expected input and output formats is provided within the package.
> 
> 
