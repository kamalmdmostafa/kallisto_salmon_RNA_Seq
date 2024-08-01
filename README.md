# RNA-Seq Quantification Pipelines

This repository contains two flexible RNA-Seq quantification pipelines using Kallisto and Salmon, as well as comparison scripts to analyze outputs from both tools.

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
   - [Installing Conda and Mamba](#installing-conda-and-mamba)
   - [Setting up the Environment](#setting-up-the-environment)
4. [Usage](#usage)
   - [Kallisto Pipeline](#kallisto-pipeline)
   - [Salmon Pipeline](#salmon-pipeline)
   - [Comparison Scripts](#comparison-scripts)
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

Additionally, comparison scripts are provided to analyze and visualize the differences between Kallisto and Salmon results.

## Requirements

- Bash (version 4.0 or later)
- Conda or Mamba (for easy installation of Kallisto and Salmon)
- Python 3.6+ (for the comparison scripts)
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
   chmod +x complete-kallisto-pipeline.sh complete-salmon-pipeline.sh
   ```

Now you have an environment with Kallisto, Salmon, and all the necessary Python libraries installed.

## Usage

Remember to activate your environment before running the pipelines:

```
conda activate rnaseq_pipelines
```

### Kallisto Pipeline

```
./complete-kallisto-pipeline.sh -cds <CDS file> -k <k-mer number> -threads <number of threads> [-o <output directory>] [-r <read file directory>] [-b <bootstrap samples>]
```

Options:
- `-cds`: Path to the CDS FASTA file
- `-k`: K-mer size for indexing
- `-threads`: Number of threads to use
- `-o`: (Optional) Output directory (default: './kallisto')
- `-r`: (Optional) Directory containing read files (default: current directory)
- `-b`: (Optional) Number of bootstrap samples (default: 100)

Example:
```
./complete-kallisto-pipeline.sh -cds path/to/cds.fasta -k 31 -threads 8 -o /path/to/custom/output -r /path/to/read/files -b 200
```

### Salmon Pipeline

```
./complete-salmon-pipeline.sh -cds <CDS file> -k <k-mer number> -threads <number of threads> [-genome <genome file>] [-vb] [-o <output directory>] [-r <read file directory>] [-b <bootstrap samples>]
```

Options:
- `-cds`: Path to the CDS FASTA file
- `-k`: K-mer size for indexing
- `-threads`: Number of threads to use
- `-genome`: (Optional) Path to the genome FASTA file for decoy-aware indexing
- `-vb`: (Optional) Use Variational Bayesian Optimization in quantification
- `-o`: (Optional) Output directory (default: './salmon')
- `-r`: (Optional) Directory containing read files (default: current directory)
- `-b`: (Optional) Number of bootstrap samples (default: 200)

Example:
```
./complete-salmon-pipeline.sh -cds path/to/cds.fasta -k 31 -threads 8 -genome path/to/genome.fasta -vb -o /path/to/custom/output -r /path/to/read/files -b 300
```

### Comparison Scripts

After running both the Kallisto and Salmon pipelines, you can use the comparison scripts to analyze the results:

1. Basic comparison script:
```
python3 compare_quantification_tools_kallisto_and_salmon.py
```

2. Alternative comparison script with additional analyses:
```
python3 compare_quantification_tools_kallisto_and_salmon_alternative.py
```

Both scripts will prompt you to enter:
1. The path to the first tool's results directory
2. The path to the second tool's results directory
3. The path to save the comparison results

The scripts will then generate various comparisons and visualizations. The alternative version produces some extra comparisons and may take longer to run, depending on the number of quantification files in your directories.

## Output

Both pipelines will create a directory (specified by the `-o` option or defaulting to `kallisto` or `salmon`) containing:

1. Cleaned CDS FASTA file
2. Index files
3. Quantification results for each sample
4. QC report (`qc_report.txt`)
5. TPM matrix (`kallisto_tpm_matrix.tsv` or `salmon_tpm_matrix.tsv`)
6. CPM matrix (`kallisto_cpm_matrix.tsv` or `salmon_cpm_matrix.tsv`)

The comparison scripts will generate:

1. An Excel file with detailed comparison statistics
2. Various plots including scatter plots, density plots, Bland-Altman plots, and MA plots
3. Correlation heatmaps
4. Transcript detection comparison plots
5. Mapping rate comparison plots
6. Violin plots for TPM, Count, and EffLength distributions
7. Cumulative TPM distribution plot

## Troubleshooting

- Ensure you have activated the environment (`conda activate rnaseq_pipelines`) before running the pipelines or comparison scripts.
- If you encounter issues with Conda, try using Mamba for faster and more reliable package installation.
- Check that input files are in the correct format (FASTA for CDS/genome, gzipped FASTQ for reads).
- For any errors, check the log files in the output directory.
- If you encounter permission issues, make sure the scripts are executable (`chmod +x script_name.sh`).
- For issues with the comparison scripts, ensure all required Python libraries are installed and up-to-date in your environment.
- If the scripts can't find your read files, make sure you're using the `-r` option to specify the correct directory.
- If you're having memory issues during the comparison, try running the script on a subset of your data first.

## Contact

For any questions or issues, please open an issue on this GitHub repository.
