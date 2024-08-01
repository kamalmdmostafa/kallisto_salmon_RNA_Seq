#!/bin/bash

set -e

# Function to display usage
usage() {
    echo "Usage: $0 -cds <CDS file> -k <k-mer number> -threads <number of threads> [-o <output directory>] [-r <read file directory>] [-b <bootstrap samples>]"
    echo "  -cds      : Path to the CDS FASTA file"
    echo "  -k        : K-mer size for indexing"
    echo "  -threads  : Number of threads to use"
    echo "  -o        : Output directory (optional, default: './kallisto')"
    echo "  -r        : Directory containing read files (optional, default: current directory)"
    echo "  -b        : Number of bootstrap samples (optional, default: 100)"
    exit 1
}

# Parse command-line arguments
OUTPUT_DIR="./kallisto"
READ_DIR="."
BOOTSTRAP_SAMPLES=100
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
        -o)
        OUTPUT_DIR="$2"
        shift 2
        ;;
        -r)
        READ_DIR="$2"
        shift 2
        ;;
        -b)
        BOOTSTRAP_SAMPLES="$2"
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

# Check if read directory exists
if [ ! -d "$READ_DIR" ]; then
    echo "Error: Read directory '$READ_DIR' does not exist"
    exit 1
fi

# Set up directory structure
KALLISTO_DIR="$OUTPUT_DIR"
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
echo "Using k-mer size: $K_VALUE"

echo "Command: kallisto index -i $INDEX_PREFIX -k $K_VALUE $CLEANED_CDS"
kallisto index -i "$INDEX_PREFIX" -k "$K_VALUE" "$CLEANED_CDS"

if [ $? -ne 0 ]; then
    echo "Error: Kallisto index creation failed with exit code $?"
    exit 1
fi

if [ -f "$INDEX_PREFIX" ]; then
    echo "Index created successfully: $INDEX_PREFIX"
else
    echo "Error: Index creation failed"
    exit 1
fi

# Step 3: Quantification
echo "Starting quantification..."
mkdir -p "$QUANT_PREFIX"

# Function to find paired read files
find_paired_reads() {
    local r1_files=($(find "$READ_DIR" -maxdepth 1 -name "*R1*.fastq.gz" -o -name "*_1.fastq.gz"))
    local r2_files=($(find "$READ_DIR" -maxdepth 1 -name "*R2*.fastq.gz" -o -name "*_2.fastq.gz"))
    
    if [ ${#r1_files[@]} -eq 0 ] || [ ${#r2_files[@]} -eq 0 ]; then
        echo "Error: No paired read files found in the directory: $READ_DIR" >&2
        return 1
    fi
    
    for r1_file in "${r1_files[@]}"; do
        local r2_file="${r1_file/R1/R2}"
        r2_file="${r2_file/_1/_2}"
        
        if [ -f "$r2_file" ]; then
            echo "$r1_file $r2_file"
        else
            echo "Warning: No matching R2 file found for $r1_file" >&2
        fi
    done
}

# Process each pair of read files
while read -r r1_file r2_file; do
    if [ -z "$r1_file" ] || [ -z "$r2_file" ]; then
        continue
    fi
    
    SAMPLE_NAME=$(basename "$r1_file" | sed -E 's/_R1.*//;s/_1\.fastq\.gz//')
    echo "Processing sample: $SAMPLE_NAME"
    echo "R1 file: $r1_file"
    echo "R2 file: $r2_file"
    
    # Retry mechanism
    MAX_RETRIES=3
    for ((i=1; i<=MAX_RETRIES; i++)); do
        echo "Attempt $i of $MAX_RETRIES"
        
        kallisto quant -i "$INDEX_PREFIX" \
                       -o "${QUANT_PREFIX}/${SAMPLE_NAME}" \
                       -t "$THREADS" \
                       --bootstrap-samples="$BOOTSTRAP_SAMPLES" \
                       --plaintext \
                       "$r1_file" "$r2_file" 2>&1 | tee "${QUANT_PREFIX}/${SAMPLE_NAME}_kallisto.log"
        
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
done < <(find_paired_reads)

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
    
    # Find all abundance.tsv files
    ABUNDANCE_FILES=($(find "$QUANT_PREFIX" -name "abundance.tsv"))
    
    if [ ${#ABUNDANCE_FILES[@]} -eq 0 ]; then
        echo "Error: No abundance.tsv files found in $QUANT_PREFIX"
        return 1
    fi
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "${ABUNDANCE_FILES[0]}" > "$TPM_MATRIX"
    
    # Extract and append TPM values for each sample
    for ABUNDANCE_FILE in "${ABUNDANCE_FILES[@]}"; do
        SAMPLE_NAME=$(basename $(dirname "$ABUNDANCE_FILE"))
        echo "Processing $SAMPLE_NAME"
        if [ -f "$ABUNDANCE_FILE" ]; then
            awk -v OFS='\t' -v sample="$SAMPLE_NAME" 'NR>1 {print $5}' "$ABUNDANCE_FILE" | \
            paste "$TPM_MATRIX" - > "${TPM_MATRIX}.tmp" && mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
        else
            echo "Warning: abundance.tsv not found for $SAMPLE_NAME"
        fi
    done
    
    # Add header with sample names
    SAMPLE_NAMES=$(for f in "${ABUNDANCE_FILES[@]}"; do basename $(dirname "$f"); done | tr '\n' '\t' | sed 's/\t$//')
    (echo -e "gene_id\t$SAMPLE_NAMES" && cat "$TPM_MATRIX") > "${TPM_MATRIX}.tmp" && \
    mv "${TPM_MATRIX}.tmp" "$TPM_MATRIX"
    
    echo "TPM matrix created: $TPM_MATRIX"
}

# Function to extract and calculate CPM values
extract_cpm() {
    echo "Extracting and calculating CPM values..."
    CPM_MATRIX="${KALLISTO_DIR}/kallisto_cpm_matrix.tsv"
    
    # Find all abundance.tsv files
    ABUNDANCE_FILES=($(find "$QUANT_PREFIX" -name "abundance.tsv"))
    
    if [ ${#ABUNDANCE_FILES[@]} -eq 0 ]; then
        echo "Error: No abundance.tsv files found in $QUANT_PREFIX"
        return 1
    fi
    
    # Create header with gene names
    awk 'NR>1 {print $1}' "${ABUNDANCE_FILES[0]}" > "$CPM_MATRIX"
    
    # Extract est_counts, calculate CPM, and append for each sample
    for ABUNDANCE_FILE in "${ABUNDANCE_FILES[@]}"; do
        SAMPLE_NAME=$(basename $(dirname "$ABUNDANCE_FILE"))
        echo "Processing $SAMPLE_NAME"
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
    SAMPLE_NAMES=$(for f in "${ABUNDANCE_FILES[@]}"; do basename $(dirname "$f"); done | tr '\n' '\t' | sed 's/\t$//')
    (echo -e "gene_id\t$SAMPLE_NAMES" && cat "$CPM_MATRIX") > "${CPM_MATRIX}.tmp" && \
    mv "${CPM_MATRIX}.tmp" "$CPM_MATRIX"
    
    echo "CPM matrix created: $CPM_MATRIX"
}

# Extract TPM and CPM matrices
extract_tpm
extract_cpm

echo "TPM and CPM matrices have been generated."
echo "Pipeline completed successfully. Results are in $KALLISTO_DIR"
