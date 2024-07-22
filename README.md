# RNA-Seq Quantification Pipelines

This repository contains two flexible RNA-Seq quantification pipelines using Kallisto and Salmon, as well as a comparison script to analyze the results from both tools.

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
   - [Installing Conda and Mamba](#installing-conda-and-mamba)
   - [Setting up the Environment](#setting-up-the-environment)
4. [Usage](#usage)
   - [Kallisto Pipeline](#kallisto-pipeline)
   - [Salmon Pipeline](#salmon-pipeline)
   - [Comparison Script](#comparison-script)
5. [Output](#output)
6. [Troubleshooting](#troubleshooting)
7. [Contact](#contact)

## Overview

These pipelines automate the process of RNA-Seq quantification using either Kallisto or Salmon. They include steps for:

1. Cleaning and checking input CDS FASTA files
2. Building indices
3. Quantifying gene expression
4. Generating quality control reports
5. Creating TPM (Transcripts Per Million) and CPM (Counts Per Million) matrices

Additionally, a comparison script is provided to analyze and visualize the differences between Kallisto and Salmon results.

## Requirements

- Bash (version 4.0 or later)
- Conda or Mamba (for easy installation of Kallisto and Salmon)
- Python 3.6+ (for the comparison script)
- Python libraries: pandas, matplotlib, seaborn, scipy, scikit-learn, numpy, biopython

## Installation

### Installing Conda and Mamba

If you don't have Conda installed:

1. Download and install Miniconda from the [official website](https://docs.conda.io/en/latest/miniconda.html).

2. After installation, initialize Conda for your shell:
   ```
   conda init
   ```

3. Close and reopen your terminal for the changes to take effect.

To install Mamba (a faster alternative to Conda):

1. Install Mamba in your base Conda environment:
   ```
   conda install mamba -n base -c conda-forge
   ```

### Setting up the Environment

Choose either Conda or Mamba for the following steps. Mamba is generally faster.

1. Create a new environment for the RNA-Seq pipelines:

   Using Conda:
   ```
   conda create -n rnaseq_pipelines python=3.9
   ```

   Using Mamba:
   ```
   mamba create -n rnaseq_pipelines python=3.9
   ```

2. Activate the new environment:
   ```
   conda activate rnaseq_pipelines
   ```

3. Install Kallisto and Salmon in the environment:

   Using Conda:
   ```
   conda install -c bioconda kallisto salmon
   ```

   Using Mamba:
   ```
   mamba install -c bioconda kallisto salmon
   ```

4. Install the required Python libraries:

   Using Conda:
   ```
   conda install pandas matplotlib seaborn scipy scikit-learn numpy biopython
   ```

   Using Mamba:
   ```
   mamba install pandas matplotlib seaborn scipy scikit-learn numpy biopython
   ```

5. Clone this repository:
   ```
   git clone https://github.com/kamalmdmostafa/kallisto_salmon_RNA_Seq.git
   cd kallisto_salmon_RNA_Seq
   ```

6. Make the pipeline scripts executable:
   ```
   chmod +x flexible_kallisto_pipeline.sh flexible_salmon_pipeline_qc.sh
   ```

Now you have an environment with Kallisto, Salmon, and all the necessary Python libraries installed.

## Usage

Remember to activate your environment before running the pipelines:

```
conda activate rnaseq_pipelines
```

### Kallisto Pipeline

```
./flexible_kallisto_pipeline.sh -cds <CDS file> -k <k-mer number> -threads <number of threads>
```

Options:
- `-cds`: Path to the CDS FASTA file
- `-k`: K-mer size for indexing
- `-threads`: Number of threads to use

Example:
```
./flexible_kallisto_pipeline.sh -cds path/to/cds.fasta -k 31 -threads 8
```

### Salmon Pipeline

```
./flexible_salmon_pipeline_qc.sh -cds <CDS file> -k <k-mer number> -threads <number of threads> [-genome <genome file>] [-vb]
```

Options:
- `-cds`: Path to the CDS FASTA file
- `-k`: K-mer size for indexing
- `-threads`: Number of threads to use
- `-genome`: (Optional) Path to the genome FASTA file for decoy-aware indexing
- `-vb`: (Optional) Use Variational Bayesian Optimization in quantification

Example:
```
./flexible_salmon_pipeline_qc.sh -cds path/to/cds.fasta -k 31 -threads 8 -genome path/to/genome.fasta -vb
```

### Comparison Script

After running both Kallisto and Salmon pipelines, you can use the comparison script to analyze the results:

```
python compare_quantification_tools_kallisto_and_salmon.py
```

The script will prompt you to enter:
1. The path to the Kallisto results directory
2. The path to the Salmon results directory
3. The path to save the comparison results

The script will then generate various comparisons and visualizations.

## Output

Both pipelines will create a directory (`kallisto` or `salmon`) containing:

1. Cleaned CDS FASTA file
2. Index files
3. Quantification results for each sample
4. QC report (`qc_report.txt`)
5. TPM matrix (`kallisto_tpm_matrix.tsv` or `salmon_tpm_matrix.tsv`)
6. CPM matrix (`kallisto_cpm_matrix.tsv` or `salmon_cpm_matrix.tsv`)

The comparison script will generate:

1. An Excel file with detailed comparison statistics
2. Various plots including scatter plots, density plots, and Bland-Altman plots
3. Correlation heatmaps
4. Transcript detection comparison plots
5. Mapping rate comparison plots

## Troubleshooting

- Ensure you have activated the environment (`conda activate rnaseq_pipelines`) before running the pipelines or comparison script.
- If you encounter issues with Conda, try using Mamba for faster and more reliable package installation.
- Check that input files are in the correct format (FASTA for CDS/genome, gzipped FASTQ for reads).
- For any errors, check the log files in the output directory.
- If you encounter permission issues, make sure the scripts are executable (`chmod +x script_name.sh`).
- For issues with the comparison script, ensure all required Python libraries are installed and up-to-date in your environment.

## Contact

For any questions or issues, please open an issue on this GitHub repository
