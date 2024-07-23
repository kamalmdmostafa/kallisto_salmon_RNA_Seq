#!/bin/bash

set -e

# Function to display usage
usage() {
    echo "Usage: $0 -cds <CDS file> -k <k-mer number> -threads <number of threads>"
    echo "  -cds      : Path to the CDS FASTA file"
    echo "  -k        : K-mer size for indexing"
    echo "  -threads  : Number of threads to use"
    exit 1
}

# Parse command-line arguments
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
KALLISTO_DIR="kallisto"
mkdir -p "$KALLISTO_DIR"

# Set variables
INDEX_PREFIX="${KALLISTO_DIR}/kallisto_index_k${K_VALUE}"
QUANT_PREFIX="${KALLISTO_DIR}/kallisto_quant_k${K_VALUE}"
CLEANED_CDS="${KALLISTO_DIR}/cleaned_$(basename "$CDS_FILE")"
QC_REPORT="${KALLISTO_DIR}/qc_report.txt"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to determine the number of threads to use
get_thread_count() {
    local total_threads=$(nproc)
    local available_threads=$((total_threads - 2))  # Leave 2 threads for system processes
    
    if [ $available_threads -lt 1 ]; then
        echo 1
    elif [ $available_threads -gt "$THREADS" ]; then
        echo "$THREADS"
    else
        echo $available_threads
    fi
}

# Check for required commands
for cmd in kallisto awk grep sed; do
    if ! command_exists "$cmd"; then
        echo "Error: $cmd is not installed or not in PATH"
        exit 1
    fi
done

echo "Kallisto version:"
kallisto version

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

# Step 2: Build Kallisto index
echo "Building Kallisto index..."
INDEX_THREADS=$(get_thread_count)
echo "Using $INDEX_THREADS threads for indexing"

kallisto index -i "$INDEX_PREFIX" -k "$K_VALUE" --threads="$INDEX_THREADS" "$CLEANED_CDS"

if [ -f "$INDEX_PREFIX" ]; then
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
    
    # Check if both R1 and R2 files exist
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "Error: Missing paired-end file for $SAMPLE_NAME"
        continue
    fi
    
    echo "Processing sample: $SAMPLE_NAME"
    
    # Retry mechanism
    MAX_RETRIES=3
    for ((i=1; i<=MAX_RETRIES; i++)); do
        echo "Attempt $i of $MAX_RETRIES"
        
        kallisto quant -i "$INDEX_PREFIX" \
                       -o "${QUANT_PREFIX}/${SAMPLE_NAME}" \
                       -t "$THREADS" \
                       --bootstrap-samples=100 \
                       "$R1_FILE" "$R2_FILE" 2>&1 | tee "${QUANT_PREFIX}/${SAMPLE_NAME}_kallisto.log"
        
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo "Quantification successful for $SAMPLE_NAME"
            break
        else
            echo "Error in quantification for $SAMPLE_NAME (Attempt $i)"
            if [ $i -eq $MAX_RETRIES ]; then
                echo "Failed to quantify $SAMPLE_NAME after $MAX_RETRIES attempts"
                echo "See log file: ${QUANT_PREFIX}/${SAMPLE_NAME}_kallisto.log"
            else
                echo "Retrying..."
                sleep 5
            fi
        fi
    done
done

# Check if any samples were processed
if [ -z "$(ls -A "$QUANT_PREFIX")" ]; then
    echo "Error: No samples were processed successfully. Check your input files and logs."
    exit 1
fi

echo "All samples have been processed. Results are in $QUANT_PREFIX"

# Step 4: Quality Check
echo "Performing quality checks..."
echo "Quality Check Report" > "$QC_REPORT"
echo "=====================" >> "$QC_REPORT"
echo "" >> "$QC_REPORT"

for SAMPLE_DIR in "${QUANT_PREFIX}"/*
do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo "Sample: $SAMPLE_NAME" >> "$QC_REPORT"
        
        if [ -f "${SAMPLE_DIR}/run_info.json" ]; then
            # Extract QC metrics from run_info.json
            TOTAL_READS=$(jq .n_processed "${SAMPLE_DIR}/run_info.json")
            MAPPED_READS=$(jq .n_pseudoaligned "${SAMPLE_DIR}/run_info.json")
            MAPPING_RATE=$(jq .p_pseudoaligned "${SAMPLE_DIR}/run_info.json")
            
            echo "  Total reads processed: $TOTAL_READS" >> "$QC_REPORT"
            echo "  Reads pseudoaligned: $MAPPED_READS" >> "$QC_REPORT"
            echo "  Pseudoalignment rate: $MAPPING_RATE" >> "$QC_REPORT"
        else
            echo "  Warning: run_info.json not found for $SAMPLE_NAME" >> "$QC_REPORT"
        fi
        
        if [ -f "${SAMPLE_DIR}/abundance.tsv" ]; then
            # Check number of quantified genes
            QUANTIFIED_GENES=$(awk 'NR>1 {count += ($5 > 0)} END {print count}' "${SAMPLE_DIR}/abundance.tsv")
            echo "  Genes quantified (TPM > 0): $QUANTIFIED_GENES" >> "$QC_REPORT"
            
            # Calculate mean TPM
            MEAN_TPM=$(awk 'NR>1 {sum += $5} END {print sum/NR}' "${SAMPLE_DIR}/abundance.tsv")
            echo "  Mean TPM: $MEAN_TPM" >> "$QC_REPORT"
        else
            echo "  Warning: abundance.tsv not found for $SAMPLE_NAME" >> "$QC_REPORT"
        fi
        
        echo "" >> "$QC_REPORT"
    fi
done

echo "QC report generated: $QC_REPORT"

# Function to extract TPM values
extract_tpm() {
    echo "Extracting TPM values..."
    TPM_MATRIX="${KALLISTO_DIR}/kallisto_tpm_matrix.tsv"
    
    # Find the first abundance.tsv file to get gene names
    FIRST_ABUNDANCE=$(find "$QUANT_PREFIX" -name "abundance.tsv" | head -n 1)
    if [ -z "$FIRST_ABUNDANCE" ]; then
        echo "Error: No abundance.tsv files found in $QUANT_PREFIX"
        return 1
    fi
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "$FIRST_ABUNDANCE" > "$TPM_MATRIX"
    
    # Extract and append TPM values for each sample
    find "$QUANT_PREFIX" -type d -mindepth 1 -maxdepth 1 | while read SAMPLE_DIR; do
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        ABUNDANCE_FILE="$SAMPLE_DIR/abundance.tsv"
        if [ -f "$ABUNDANCE_FILE" ]; then
            awk -v OFS='\t' -v sample="$SAMPLE_NAME" 'NR>1 {print $5}' "$ABUNDANCE_FILE" | \
            paste "$TPM_MATRIX" - > "${TPM_MATRIX}.tmp" && mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
        else
            echo "Warning: abundance.tsv not found for $SAMPLE_NAME"
        fi
    done
    
    # Add header with sample names
    (echo -e "gene_id\t$(find "$QUANT_PREFIX" -type d -mindepth 1 -maxdepth 1 | xargs -n 1 basename | tr '\n' '\t' | sed 's/\t$//')" && cat "$TPM_MATRIX") > "${TPM_MATRIX}.tmp" && \
    mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
    
    echo "TPM matrix created: $TPM_MATRIX"
}

# Function to extract and calculate CPM values
extract_cpm() {
    echo "Extracting and calculating CPM values..."
    CPM_MATRIX="${KALLISTO_DIR}/kallisto_cpm_matrix.tsv"
    
    # Find the first abundance.tsv file to get gene names
    FIRST_ABUNDANCE=$(find "$QUANT_PREFIX" -name "abundance.tsv" | head -n 1)
    if [ -z "$FIRST_ABUNDANCE" ]; then
        echo "Error: No abundance.tsv files found in $QUANT_PREFIX"
        return 1
    }
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "$FIRST_ABUNDANCE" > "$CPM_MATRIX"
    
    # Extract est_counts, calculate CPM, and append for each sample
    find "$QUANT_PREFIX" -type d -mindepth 1 -maxdepth 1 | while read SAMPLE_DIR; do
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        ABUNDANCE_FILE="$SAMPLE_DIR/abundance.tsv"
        if [ -f "$ABUNDANCE_FILE" ]; then
            awk -v OFS='\t' -v sample="$SAMPLE_NAME" '
            NR>1 {
                count[NR] = $4
                sum += $4
            }
            END {
                for (i in count) {
                    print (count[i] / sum) * 1e6
                }
            }' "$ABUNDANCE_FILE" | \
            paste "$CPM_MATRIX" - > "${CPM_MATRIX}.tmp" && mv "${CPM_MATRIX}.tmp" "$CPM_MATRIX"
        else
            echo "Warning: abundance.tsv not found for $SAMPLE_NAME"
        fi
    done
    
    # Add header with sample names
    (echo -e "gene_id\t$(find "$QUANT_PREFIX" -type d -mindepth 1 -maxdepth 1 | xargs -n 1 basename | tr '\n' '\t' | sed 's/\t$//')" && cat "$CPM_MATRIX") > "${CPM_MATRIX}.tmp" && \
    mv "${CPM_MATRIX}.tmp" "$CPM_MATRIX"
    
    echo "CPM matrix created: $CPM_MATRIX"
}

# Extract TPM and CPM matrices
extract_tpm
extract_cpm

echo "TPM and CPM matrices have been generated."
echo "Pipeline completed successfully."
