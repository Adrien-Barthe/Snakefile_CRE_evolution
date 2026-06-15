# Modernized Pipeline for CRE Evolution

## Project Overview
The goal of this project was to modernize an older pipeline that relied on multiple separated Bash scripts. It has now been unified into a single Snakefile. 

This pipeline automates the entire workflow: it generates genomic models, runs SLiM simulations in parallel, processes raw outputs, calculates statistics, and generates plots.

## Repository Structure
In order to keep this repository as clean and light as possible, data and R scripts are ignored by git.

* `Snakefile` : The core Snakemake pipeline containing all rules and R script recipes.
* `environment.yml` : Conda environment file containing all necessary dependencies.
* `.gitignore` : Keeps the repo clean by ignoring raw data and dynamically generated folders.

*(Note: All folders are generated automatically by the pipeline during running).*

## Prerequisites & Installation
This pipeline is designed to run on an HPC cluster (using SLURM). To guarantee reproducibility, an `environment.yml` file is provided.

First clone the repository:
```bash
git clone https://github.com/Adrien-Barthe/Snakefile_CRE_evolution.git
```

Second create and activate the Conda environment:
```bash
cd Snakefile_CRE_evolution
conda env create -f environment.yml
conda activate pipeline_reg
```
To ensure the pipeline runs smoothly on any machine or cluster, we redirect the temporary files to a local `tmp/` folder directly inside the current project directory:

```bash
# Creates a local tmp folder and tells SLiM to use it
export TMPDIR=$PWD/tmp
mkdir -p $TMPDIR
```

## Pipeline Architecture
The Snakemake workflow handles the following steps automatically:
* **create_models**: Writes the specific .txt Eidos scripts for epistatic and null models.
* **run_slim_***: Executes the simulations using xargs -P for efficient threading.
* **merge_***: Aggregates the mutation counts across all replicates.
* **create_R_scripts**: Writes all required R scripts onto 02-scripts/. This guarantees that the R script is up to date with the rest of the snakefile.
* **compute_sum_stats**: Uses a dynamic batching system to parse SLiM outputs and compute pi, perSite, and haps stats safely.
* **plot_ rules**: Launches the R scripts to perform Wilcoxon tests, extract p-values, and draw ggpubr boxplots.

## How to Run the Pipeline

The pipeline has two modes, selected via `--config plot_only=...`:

* **Full mode (default)**: runs SLiM simulations for the parameters defined at the top of the Snakefile (if needed), computes summary statistics, and produces plots from both the new simulations and any previously computed data.
```bash
Launch the pipeline on the cluster ( adjust jobs parameters accordingly to cluster capacity, 5 jobs means that at maximum you take 100 cores and 160 GB ram)
nohup snakemake --executor slurm --jobs 5 &
```

* **Plot-only mode**: never launches simulations or stat computations. It only regenerates PDFs from data that already exists in `05-results/`.
```bash
nohup snakemake --config plot_only=True --executor slurm --jobs 5 &
```

Since plot-only mode doesn't need any cluster resources, it can also be run locally on any machine as long as the conda environment is installed and activated, and `05-results/is present:
```bash
snakemake --config plot_only=True --cores 1
```


## Outputs & Results
Once finished, all your results will be neatly organized inside the 05-results/ directory:
* **pdfs_sm2_analysis/**: SFS, haplotypes, and fixed stats plots.
* **pdfs_sm3_analysis/**: Comparisons based on sm3 coefficients.
* **pdfs_h_analysis/**: Dominance effect plots.
* **pdfs_rec_analysis/**: Recombination and physical distance plots.
* **p_values/**: Exact p-values (.txt) from Wilcoxon tests corresponding to the plots.
* **sumStats/**: Means, medians, and Fisher Exact Test results.

## Configuration & Smart Execution
At the top of the Snakefile, you will find the Simulation Parameters. You only need to input the parameters you want to simulate now.

```python
# Example: What SLiM will run right now
DOMINANCES = ["0.1", "0.7"] 
SELECTIONS = ["-0.001", "-0.01", "-0.0425", "-0.1", "-0.3"]
SELECTIONS_M3 = ["-0.0001"]

```

* **DOMINANCES**: List of dominance coefficients (h) to test.
* **SELECTIONS**: Selection coefficients for one type of coding mutations (sM2).
* **SELECTIONS_M3**: Selection coefficients for regulatory mutations (sM3).
* **REPLICATS**: The number (-1) of replicates per model (ex: range(1, 1001)=1000 simulations).
* **GENERATION**: The number of generations before outputting data (default: 10000).
* **CRE_L**: Cis-Regulatory Element length in kb.
* **PD_L**: Physical distance between genomic elements in kb.

A second block, the Plot Locking System, defines the full ranges a plot needs before it is unlocked (before Snakemake will request it):

```python
TARGET_H = ["0.1", "0.25", "0.5", "0.7"]
TARGET_REC = ["0.5", "1", "5", "10", "500"]
TARGET_SM3 = ["-0", "-0.0001", "-0.001"]
TARGET_S = ["-0.001", "-0.01", "-0.0425", "-0.1", "-0.3"]
```

* **TARGET_H / TARGET_REC / TARGET_SM3 / TARGET_S**: the full set of h / recombination distance / sM3 / sM2 values that a given comparison plot spans across. A plot is only added to the run's targets once all the required combinations exist (from current simulations, past simulations, or both).

**Smart Plotting Auto-Detection**: You do not need to rerun old simulations to plot them together with new ones. At startup, the Snakefile scans `05-results/` to find which parameter combinations already exist. Combinations from `DOMINANCES`/`SELECTIONS`/`SELECTIONS_M3`/`PD_L`/`CRE_L` that are already present are skipped entirely (no simulation, stats, or model files are regenerated for them), only new combinations are simulated. The plot-locking system then considers both the new and the pre-existing combinations together, so comparison plots can mix old and new data.
