#!/bin/bash

set -e

# Function to display usage
usage() {
    echo "Usage: $0 -cds <CDS file> -k <k-mer number> -threads <number of threads> [-genome <genome file>] [-vb] [-o <output directory>]"
    echo "  -cds      : Path to the CDS FASTA file"
    echo "  -k        : K-mer size for indexing"
    echo "  -threads  : Number of threads to use"
    echo "  -genome   : (Optional) Path to the genome FASTA file for decoy-aware indexing"
    echo "  -vb       : (Optional) Use Variational Bayesian Optimization in quantification"
    echo "  -o        : (Optional) Output directory (default: './salmon')"
    exit 1
}

# Parse command-line arguments
USE_VBOPT=false
OUTPUT_DIR="./salmon"
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -cds)
        CDS_FILE="$2"
        shift 2
        ;;
        -k)
        K_VALUE="$2"
        shift 2
        ;;
        -threads)
        THREADS="$2"
        shift 2
        ;;
        -genome)
        GENOME_FILE="$2"
        shift 2
        ;;
        -vb)
        USE_VBOPT=true
        shift
        ;;
        -o)
        OUTPUT_DIR="$2"
        shift 2
        ;;
        *)
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

# Check if required arguments are provided
if [ -z "$CDS_FILE" ] || [ -z "$K_VALUE" ] || [ -z "$THREADS" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Set up directory structure
SALMON_DIR="$OUTPUT_DIR"
mkdir -p "$SALMON_DIR"

# Set variables based on whether a genome file is provided
if [ -n "$GENOME_FILE" ]; then
    INDEX_PREFIX="${SALMON_DIR}/salmon_index_decoy_k${K_VALUE}"
    QUANT_PREFIX="${SALMON_DIR}/salmon_quant_decoy_k${K_VALUE}"
else
    INDEX_PREFIX="${SALMON_DIR}/salmon_index_cds_only_k${K_VALUE}"
    QUANT_PREFIX="${SALMON_DIR}/salmon_quant_cds_only_k${K_VALUE}"
fi

CLEANED_CDS="${SALMON_DIR}/cleaned_$(basename "$CDS_FILE")"
COMBINED_REF="${SALMON_DIR}/gentrome.fa"
DECOY_LIST="${SALMON_DIR}/decoys.txt"
QC_REPORT="${SALMON_DIR}/qc_report.txt"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required commands
for cmd in salmon awk grep sed; do
    if ! command_exists "$cmd"; then
        echo "Error: $cmd is not installed or not in PATH"
        exit 1
    fi
done

echo "Salmon version:"
salmon --version

# Step 1: Check and clean CDS FASTA
echo "Checking and cleaning CDS FASTA..."
awk '
BEGIN {OFS="\n"}
/^>/ {
    if (NR>1 && len>=20) 
        print header, seq
    header=$0
    seq=""
    len=0
    next
}
{
    gsub(/[^ATGCN]/, "N")
    seq=seq $0
    len+=length($0)
}
END {
    if (len>=20)
        print header, seq
}' "$CDS_FILE" > "$CLEANED_CDS"

echo "Cleaned CDS file created: $CLEANED_CDS"

# Step 2: Build Salmon index
echo "Building Salmon index..."
if [ -n "$GENOME_FILE" ]; then
    echo "Creating decoy-aware reference..."
    cat "$CLEANED_CDS" "$GENOME_FILE" > "$COMBINED_REF"
    grep "^>" "$GENOME_FILE" | cut -d " " -f 1 | sed 's/>//' > "$DECOY_LIST"
    
    salmon index -t "$COMBINED_REF" \
                 -d "$DECOY_LIST" \
                 -i "$INDEX_PREFIX" \
                 -k "$K_VALUE" \
                 --keepDuplicates \
                 -p "$THREADS"
else
    salmon index -t "$CLEANED_CDS" \
                 -i "$INDEX_PREFIX" \
                 -k "$K_VALUE" \
                 --keepDuplicates \
                 -p "$THREADS"
fi

if [ -d "$INDEX_PREFIX" ]; then
    echo "Index created successfully: $INDEX_PREFIX"
else
    echo "Error: Index creation failed"
    exit 1
fi

# Step 3: Quantification
echo "Starting quantification..."
mkdir -p "$QUANT_PREFIX"

for R1_FILE in *_R1.fastq.gz
do
    R2_FILE="${R1_FILE/_R1/_R2}"
    SAMPLE_NAME=$(basename "$R1_FILE" _R1.fastq.gz)
    
    echo "Processing sample: $SAMPLE_NAME"
    
    VBOPT_FLAG=""
    if [ "$USE_VBOPT" = true ]; then
        VBOPT_FLAG="--useVBOpt"
    fi
    
    salmon quant -i "$INDEX_PREFIX" \
                 -l A \
                 -1 "$R1_FILE" \
                 -2 "$R2_FILE" \
                 -p "$THREADS" \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 --rangeFactorizationBins 4 \
                 --numBootstraps 200 \
                 $VBOPT_FLAG \
                 -o "${QUANT_PREFIX}/${SAMPLE_NAME}"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All samples have been processed. Results are in $QUANT_PREFIX"

# Step 4: Quality Check
echo "Performing quality checks..."
echo "Quality Check Report" > "$QC_REPORT"
echo "=====================" >> "$QC_REPORT"
echo "" >> "$QC_REPORT"

for SAMPLE_DIR in "${QUANT_PREFIX}"/*
do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    echo "Sample: $SAMPLE_NAME" >> "$QC_REPORT"
    
    # Extract QC metrics from logs
    MAPPED_READS=$(grep "Mapping rate" "${SAMPLE_DIR}/logs/salmon_quant.log" | awk '{print $8}')
    MAPPING_RATE=$(grep "Mapping rate" "${SAMPLE_DIR}/logs/salmon_quant.log" | awk '{print $10}')
    
    echo "  Reads mapped: $MAPPED_READS" >> "$QC_REPORT"
    echo "  Mapping rate: $MAPPING_RATE" >> "$QC_REPORT"
    
    # Check number of quantified genes
    QUANTIFIED_GENES=$(awk 'NR>1 && $4 > 0 {count++} END {print count}' "${SAMPLE_DIR}/quant.sf")
    echo "  Genes quantified (TPM > 0): $QUANTIFIED_GENES" >> "$QC_REPORT"
    
    # Calculate mean TPM
    MEAN_TPM=$(awk 'NR>1 {sum += $4} END {print sum/(NR-1)}' "${SAMPLE_DIR}/quant.sf")
    echo "  Mean TPM: $MEAN_TPM" >> "$QC_REPORT"
    
    echo "" >> "$QC_REPORT"
done

echo "QC report generated: $QC_REPORT"

# Function to extract TPM values
extract_tpm() {
    echo "Extracting TPM values..."
    TPM_MATRIX="${SALMON_DIR}/salmon_tpm_matrix.tsv"
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "${QUANT_PREFIX}/$(ls "${QUANT_PREFIX}" | head -n1)/quant.sf" > "$TPM_MATRIX"
    
    # Extract and append TPM values for each sample
    for SAMPLE_DIR in "${QUANT_PREFIX}"/*; do
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        awk -v OFS='\t' -v sample="$SAMPLE_NAME" 'NR>1 {print $4}' "${SAMPLE_DIR}/quant.sf" | \
        paste "$TPM_MATRIX" - > "${TPM_MATRIX}.tmp" && mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
    done
    
    # Add header with sample names
    (echo -e "gene_id\t$(ls "${QUANT_PREFIX}" | tr '\n' '\t' | sed 's/\t$//')" && cat "$TPM_MATRIX") > "${TPM_MATRIX}.tmp" && \
    mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
    
    echo "TPM matrix created: $TPM_MATRIX"
}

# Function to extract and calculate CPM values
extract_cpm() {
    echo "Extracting and calculating CPM values..."
    CPM_MATRIX="${SALMON_DIR}/salmon_cpm_matrix.tsv"
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "${QUANT_PREFIX}/$(ls "${QUANT_PREFIX}" | head -n1)/quant.sf" > "$CPM_MATRIX"
    
    # Extract NumReads, calculate CPM, and append for each sample
    for SAMPLE_DIR in "${QUANT_PREFIX}"/*; do
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        awk -v OFS='\t' -v sample="$SAMPLE_NAME" '
        NR>1 {
            count[NR] = $5
            sum += $5
        }
        END {
            for (i in count) {
                print (count[i] / sum) * 1e6
            }
        }' "${SAMPLE_DIR}/quant.sf" | \
        paste "$CPM_MATRIX" - > "${CPM_MATRIX}.tmp" && mv "${CPM_MATRIX}.tmp" "$CPM_MATRIX"
    done
    
    # Add header with sample names
    (echo -e "gene_id\t$(ls "${QUANT_PREFIX}" | tr '\n' '\t' | sed 's/\t$//')" && cat "$CPM_MATRIX") > "${CPM_MATRIX}.tmp" && \
    mv "${CPM_MATRIX}.tmp" "$CPM_MATRIX"
    
    echo "CPM matrix created: $CPM_MATRIX"
}

# Extract TPM and CPM matrices
extract_tpm
extract_cpm

echo "TPM and CPM matrices have been generated."
echo "Pipeline completed successfully. Results are in $SALMON_DIR"
