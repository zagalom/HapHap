#!/bin/bash

# File: phasing_and_haplotype_pipeline.sh
# Purpose: Automate mapping, phasing (GATK + WhatsHap) and
# haplotype reconstruction (Python) using FASTQ/BAM files.

set -euo pipefail

# --- CRITICAL LOCALE FIX: Ensures 'printf' and 'bc' use the dot (.) as a decimal separator ---
export LC_NUMERIC="C"

# ---- CONFIGURATION ----

# Allows the user to provide input or prompts if not provided (with default values)
GENES_FASTA=${1:-}
RESULTS_DIR=${2:-}
FASTQ_DIR=${3:-}
THREADS=${4:-4}      # Default: 4 threads
GENES_LIST="${5:-}"  # Optional argument, e.g. ERG11,GSC1,UPC2

# MULTIPLICATION FACTOR FOR MAX DP FILTER (MEAN_DP * FACTOR)
# We use 2X to mitigate gene collapse and only remove extreme outliers.
DP_MULTIPLIER="2"

# ATTENTION: Set this variable to "true" to SKIP Allelic Depth (AD) filtering.
BYPASS_FILTER="false"

# Read Group (RG) variables - Inserted into the header by BWA
RG_LB="LIB1"
RG_PL="ILLUMINA"
RG_PU="L001"

HAPLOTYPE_SCRIPT="haplotype_reconstruction.py"
TEMP_VCF_SUFFIX=".unphased.vcf.gz"
FILTERED_VCF_SUFFIX=".filtered.vcf.gz" # Temporary filtered VCF (PASS)
FINAL_VCF_SUFFIX=".filtered.phased.vcf.gz"
IMBALANCE_REPORT_SUFFIX=".imbalanced_report.txt"
ALL_VARIANTS_REPORT_SUFFIX=".all_variants_filter_report.txt"

function print_usage {
    echo "Usage:"
    echo "  $0 <genes_fasta> <results_dir> <fastq_dir> [threads] [genes_list]"
    echo "    <genes_fasta>  : Reference genes FASTA file"
    echo "    <results_dir>  : Directory for output/sample subdirectories"
    echo "    <fastq_dir>    : Directory containing paired FASTQs (SAMPLE_1.fastq.gz/SAMPLE_2.fastq.gz)"
    echo "    [threads]      : Number of threads to use (default: 4)"
    echo "    [genes_list]   : Optional: comma-separated genes (e.g. ERG11,GSC1). Restricts analysis."
    exit 1
}

# Checks if mandatory arguments were provided
[ -n "$GENES_FASTA" ] && [ -n "$RESULTS_DIR" ] && [ -n "$FASTQ_DIR" ] || print_usage

# Defines the region filter for GATK/bcftools (only if GENES_LIST is provided)
REGION_FILTER=""
if [ -n "$GENES_LIST" ]; then
    echo "Warning: Analysis will be RESTRICTED to genes: $GENES_LIST."
    # Transforms 'G1,G2,G3' into the format '-L G1 -L G2 -L G3' for GATK
    REGION_FILTER=$(echo "$GENES_LIST" | tr ',' '\n' | sed 's/^/-L /' | tr '\n' ' ')
fi


# --- CLEANUP AND REPORTING FUNCTIONS ---

# Function to clean up temporary files (defined here to be used in the loop)
cleanup_temp_files() {
    local ISOLATE_DIR=$1
    local ISOLATE_NAME=$2
    # Temporary variables (need to redefine them locally)
    local TEMP_DECOMPOSED="$ISOLATE_DIR/$ISOLATE_NAME.decomposed.vcf"
    local TEMP_SNPS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.snps.vcf.gz"
    local TEMP_INDELS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.indels.vcf.gz"
    local TEMP_SNPS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.snps.filtered.vcf.gz"
    local TEMP_INDELS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.indels.filtered.vcf.gz"
    local PASS_POSITIONS_FILE="$ISOLATE_DIR/$ISOLATE_NAME.pass_positions.tmp"
    local UNPHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${TEMP_VCF_SUFFIX}"

    # Cleanup of intermediate temporary files
    rm -f "$TEMP_DECOMPOSED"
    rm -f "$TEMP_SNPS_VCF" "$TEMP_SNPS_VCF.tbi"
    rm -f "$TEMP_INDELS_VCF" "$TEMP_INDELS_VCF.tbi"
    rm -f "$TEMP_SNPS_FILTERED" "$TEMP_SNPS_FILTERED.tbi"
    rm -f "$TEMP_INDELS_FILTERED" "$TEMP_INDELS_FILTERED.tbi"
    rm -f "$PASS_POSITIONS_FILE" # Cleanup of the temporary file
    rm -f ${UNPHASED_VCF}
    rm -f ${UNPHASED_VCF%.*}.tbi
}
# -------------------------------------


# ---- PHASE 0: Mapping FASTQs and generating BAMs ----

echo "### STEP 0: Mapping reads and generating sorted BAMs..."
echo "    -> Using $THREADS threads."

if [ ! -f "${GENES_FASTA}.bwt" ]; then
    echo "    -> Indexing $GENES_FASTA for BWA..."
    bwa index "$GENES_FASTA"
fi

mkdir -p "$RESULTS_DIR"

for fq1 in "$FASTQ_DIR"/*_1*.fastq*; do
    [ -e "$fq1" ] || continue
    SAMPLE=$(basename "$fq1" | sed 's/_1.*//')
    fq2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    [ -f "$fq2" ] || fq2="$FASTQ_DIR/${SAMPLE}_2.fastq"
    if [ ! -f "$fq2" ]; then
        echo "    -> Warning: R2 file (pair) not found for $fq1. Skipping sample $SAMPLE."
        continue
    fi
    smdir="$RESULTS_DIR/$SAMPLE"
    mkdir -p "$smdir"
    outbam="$smdir/${SAMPLE}.sorted.bam"
    if [ -s "$outbam" ]; then
        echo "    [$SAMPLE] Sorted BAM exists. Skipping alignment."
        continue
    fi
    echo "    [$SAMPLE] Aligning and sorting with BWA/SAMtools (@RG Header)..."
    # Adds the appropriate @RG header directly in the bwa mem command
    bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}" "$GENES_FASTA" "$fq1" "$fq2" |
        samtools view -bh -@ "$THREADS" - |
        samtools sort -@ "$THREADS" -o "$outbam"
    samtools index "$outbam"
    echo "    [$SAMPLE] BAM Complete: $outbam"
done

# ---- STEP 1: Prepare reference indices for variant calling ----

echo "### STEP 1: Preparing reference indices for GATK/BCFTOOLS..."

# 1.1. GATK requires a .dict index for the FASTA reference.
if [ ! -f "${GENES_FASTA%.*}.dict" ]; then
    echo "    -> Creating .dict index (Sequence Dictionary) for GATK..."
    gatk CreateSequenceDictionary -R ${GENES_FASTA}
fi

# 1.2. Samtools/GATK also requires the .fai index (FASTA Index)
if [ ! -f "${GENES_FASTA}.fai" ]; then
    echo "    -> Creating .fai index (FASTA Index) for Samtools/pysam..."
    samtools faidx ${GENES_FASTA}
fi


# --- LOOP PREPARATION ---
ALL_ISOLATES=$(find "$RESULTS_DIR" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep -v '^\.$' | grep -v '^..$' | grep -v "^\.$")
NUM_ISOLATES=$(echo "$ALL_ISOLATES" | wc -l)

if [ "$NUM_ISOLATES" -eq 0 ]; then
    echo "Error: No isolate subdirectories (BAMs) found in $RESULTS_DIR."
    exit 1
fi

echo "Processing $NUM_ISOLATES isolates found in $RESULTS_DIR..."
CURRENT_COUNT=0

# ---- PHASE 2: Per-sample variant calling, phasing, and haplotype reconstruction (GATK -> WhatsHap -> Python) ----

echo "### STEP 2: Calling, Filtering, and Phasing (GATK, WhatsHap, Python)"
echo "------------------------------------------------------"

for ISOLATE_NAME in $ALL_ISOLATES; do
    CURRENT_COUNT=$((CURRENT_COUNT + 1))
    ISOLATE_DIR="$RESULTS_DIR/$ISOLATE_NAME"
    INPUT_BAM_TO_USE="$ISOLATE_DIR/$ISOLATE_NAME.sorted.bam" # BAM created in STEP 0
    UNPHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${TEMP_VCF_SUFFIX}"
    FILTERED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${FILTERED_VCF_SUFFIX}" # Temporary filtered VCF
    PHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${FINAL_VCF_SUFFIX}"
    IMBALANCE_REPORT="$ISOLATE_DIR/$ISOLATE_NAME${IMBALANCE_REPORT_SUFFIX}"
    ALL_VARIANTS_REPORT="$ISOLATE_DIR/$ISOLATE_NAME${ALL_VARIANTS_REPORT_SUFFIX}"

    # Temporary file for PASS positions for reconciliation
    PASS_POSITIONS_FILE="$ISOLATE_DIR/$ISOLATE_NAME.pass_positions.tmp"

    # Temporary variables for filtering:
    TEMP_DECOMPOSED="$ISOLATE_DIR/$ISOLATE_NAME.decomposed.vcf" # VCF (uncompressed)
    TEMP_SNPS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.snps.vcf.gz"
    TEMP_INDELS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.indels.vcf.gz"
    TEMP_SNPS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.snps.filtered.vcf.gz"
    TEMP_INDELS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.indels.filtered.vcf.gz"


    echo "[$CURRENT_COUNT/$NUM_ISOLATES] -> Processing isolate: $ISOLATE_NAME"

    # Source BAM file check
    if [ ! -f "$INPUT_BAM_TO_USE" ]; then
        echo "    Warning: BAM file not found at $INPUT_BAM_TO_USE. Skipping $ISOLATE_NAME."
        continue
    fi
    # ------------------------------------------------------------------------


    # --- 2.1. VARIANT CALLING WITH GATK (HaplotypeCaller) ---
    echo "    -> 2.1. Calling variants with GATK (creates $TEMP_VCF_SUFFIX)..."
    if [ -n "$REGION_FILTER" ]; then
        echo "    -> CRITICAL: Using region filter: $GENES_LIST"
    fi

    # Use the folder name as --sample-name
    gatk HaplotypeCaller \
        -R ${GENES_FASTA} \
        -I ${INPUT_BAM_TO_USE} \
        -O ${UNPHASED_VCF} \
        ${REGION_FILTER} \
        --minimum-mapping-quality 30 \
        -ploidy 2 \
        --sample-name "${ISOLATE_NAME}" \
        --min-base-quality-score 20

    # Checks if GATK failed before proceeding
    if [ $? -ne 0 ]; then
        echo "    CRITICAL GATK ERROR for isolate $ISOLATE_NAME. Skipping phasing/reconstruction."
        cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
        continue
    fi

    # --- DEFINITION OF FILTERING THRESHOLDS (INCLUDES MIN QD AND NEW AF FILTER) ---
    SNPS_RPRS="-8.0"
    SNPS_BQRS="-8.0"
    SNPS_AD_ALT="10"
    SNPS_DP_TOTAL="20"
    SNPS_QD="2.0"

    INDELS_RPRS="-20.0"
    INDELS_BQRS="-8.0"
    INDELS_AD_ALT="10"
    INDELS_DP_TOTAL="20"
    INDELS_QD="10.0" # Changed from 25.0 to 10.0

    # NEW: Allelic Fraction (AF) limits for Heterozygotes (AF must be between 0.4 and 0.6)
    AF_MIN_BALANCE="0.4"
    AF_MAX_BALANCE="0.6"

    # DEFAULT VALUE for Maximum Depth (will be replaced by dynamic calculation)
    MAX_DP_FILTER="1000"

    # --- DYNAMIC CALCULATION OF MAX_DP_FILTER (ROBUST) ---
    echo "    -> Dynamically calculating MAX_DP_FILTER based on the isolate's average DP..."

    DP_VALUES=$(bcftools query -f '[%DP]\n' ${UNPHASED_VCF} -i 'FORMAT/DP>0' 2>/dev/null)
    NUM_VARIANTS=$(echo "$DP_VALUES" | wc -l)

    if [ "$NUM_VARIANTS" -gt 0 ]; then
        # Added || true to bc to avoid failure if the input is invalid/empty
        DP_MEAN_FLOAT=$(echo "$DP_VALUES" | awk 'BEGIN {sum=0} {sum+=$1} END {printf "%.3f", sum/NR}' || true)

        if [ -n "$DP_MEAN_FLOAT" ]; then
            NEW_MAX_DP_FLOAT=$(echo "scale=3; $DP_MEAN_FLOAT * $DP_MULTIPLIER" | bc)
            NEW_MAX_DP_FILTER=$(printf "%.0f\n" $NEW_MAX_DP_FLOAT)

            MIN_LIMIT=$((SNPS_DP_TOTAL * 2))

            if [ "$NEW_MAX_DP_FILTER" -gt "$MIN_LIMIT" ] 2>/dev/null; then
                MAX_DP_FILTER="$NEW_MAX_DP_FILTER"
                echo "        -> Calculated Mean DP: $(printf "%.2f" "$DP_MEAN_FLOAT")X"
                echo "        -> Maximum DP Limit (Rounded): ${MAX_DP_FILTER}X (Dynamic)"
            else
                echo "        -> Warning: Calculation resulted in a very low limit. Keeping default value: ${MAX_DP_FILTER}X."
            fi
        else
            echo "        -> Warning: Failed to calculate Mean DP. Keeping default value: ${MAX_DP_FILTER}X."
        fi
    else
        echo "        -> Warning: No valid variants (DP>0) to calculate Mean DP. Keeping default value: ${MAX_DP_FILTER}X."
    fi
    # -------------------------------------------------------------------------


    # --- DEFINITION OF FILTERING EXPRESSIONS ---

    # 0. NON_RANKSUM_FILTER (Minimum Depth and Quality Conditions)
    NON_RANKSUM_SNPS="FMT/AD[0:1] > ${SNPS_AD_ALT} & FMT/DP[0] > ${SNPS_DP_TOTAL} & FMT/DP[0] < ${MAX_DP_FILTER} & INFO/QD > ${SNPS_QD}"
    NON_RANKSUM_INDELS="FMT/AD[0:1] > ${INDELS_AD_ALT} & FMT/DP[0] > ${INDELS_DP_TOTAL} & FMT/DP[0] < ${MAX_DP_FILTER} & INFO/QD > ${INDELS_QD}"

    # 1. PATH A: PERFECT 1/1 (Bypass RankSum - Condition: GT=1/1 AND AD_REF=0)
    PERFECT_1_1_SNPS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & (${NON_RANKSUM_SNPS})"
    PERFECT_1_1_INDELS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & (${NON_RANKSUM_INDELS})"

    # 2. PATH B: STANDARD QUALITY CHECK (Applies all RankSums + AF Balance)
    STRICT_Q_SNPS="(ReadPosRankSum > ${SNPS_RPRS} & BaseQRankSum > ${SNPS_BQRS})"
    STRICT_Q_INDELS="(ReadPosRankSum > ${INDELS_RPRS} & BaseQRankSum > ${INDELS_BQRS})"

    STANDARD_AF_GT_CHECK="(FMT/GT=\"1/1\" | (FMT/GT=\"0/1\" & (FMT/AD[0:1] / FMT/DP[0]) >= ${AF_MIN_BALANCE} & (FMT/AD[0:1] / FMT/DP[0]) <= ${AF_MAX_BALANCE}))"

    STANDARD_PASS_SNPS="(${NON_RANKSUM_SNPS} & ${STRICT_Q_SNPS} & ${STANDARD_AF_GT_CHECK})"
    STANDARD_PASS_INDELS="(${NON_RANKSUM_INDELS} & ${STRICT_Q_INDELS} & ${STANDARD_AF_GT_CHECK})"

    # 3. FINAL FILTER: PATH A (PERFECT 1/1) OR PATH B (STANDARD PASS)
    SNPS_FINAL_FILTER="(${PERFECT_1_1_SNPS}) | (${STANDARD_PASS_SNPS})"
    INDELS_FINAL_FILTER="(${PERFECT_1_1_INDELS}) | (${STANDARD_PASS_INDELS})"

    # 4. Filter for IMBALANCE REPORT (Unbalanced 0/1 that passes Basic Quality)
    IMBALANCE_AF_CHECK_FILTER="FMT/GT=\"0/1\" & ((FMT/AD[0:1] / FMT/DP[0]) < ${AF_MIN_BALANCE} | (FMT/AD[0:1] / FMT/DP[0]) > ${AF_MAX_BALANCE})"
    # -------------------------------------------------------------------------


    # --- 2.1.5. FILTERING BY READ SUPPORT (AD) AND QUALITY (RankSums + QD) ---
    echo "    -> 2.1.5. Applying robust filtering (DP < ${MAX_DP_FILTER}X Dynamic)..."
    echo "    -> New dual path filter implemented to accept perfect 1/1 (AD_REF=0)."

    COUNT_BEFORE=$(bcftools view -H ${UNPHASED_VCF} 2>/dev/null | wc -l)
    echo "        -> Variants before filtering: ${COUNT_BEFORE} (data lines)"

    if [ "$BYPASS_FILTER" = "true" ]; then
        echo "        [DEBUG MODE] FILTERING SKIPPED (BYPASS). Simply copying..."
        # Copy and ensure BGZF and indexing
        bcftools view -O z -o "$FILTERED_VCF" ${UNPHASED_VCF}
    else
        # 2.1.5.a. Decompose multiallelic blocks and normalize
        echo "        -> 2.1.5.a. Decomposing/Normalizing variants (bcftools norm)..."
        bcftools norm -f ${GENES_FASTA} -m -any -O v -o "$TEMP_DECOMPOSED" ${UNPHASED_VCF}

        if [ $? -ne 0 ]; then
            echo "    CRITICAL BCFTOOLS NORM ERROR: Failed decomposition for isolate $ISOLATE_NAME. Skipping."
            cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
            continue
        fi

        # -------------------------------------------------------------------------
        # --- 2.1.2. GENERATING FULL FILTERING REPORT (FIRST PASS - FAIL REASONS) ---
        echo "    -> 2.1.2. Generating Full Filtering Report ($ALL_VARIANTS_REPORT_SUFFIX) with fail reasons (Status PENDING)..."

        # Filter evaluation expressions
        EVAL_SNPS_Q_FILTER="${NON_RANKSUM_SNPS}"
        EVAL_INDELS_Q_FILTER="${NON_RANKSUM_INDELS}"
        EVAL_AF_BALANCE_FILTER="FMT/GT=\"0/1\" & (FMT/AD[0:1] / FMT/DP[0]) >= ${AF_MIN_BALANCE} & (FMT/AD[0:1] / FMT/DP[0]) <= ${AF_MAX_BALANCE}"

        # Report Header
        echo -e "# Complete Variant Evaluation Report\n# ISOLATE: $ISOLATE_NAME" > "$ALL_VARIANTS_REPORT"
        echo -e "# Threshold Definitions:\n# Min DP Total: $SNPS_DP_TOTAL, Max DP Total: $MAX_DP_FILTER, Min AD Alt: $SNPS_AD_ALT, QD Min SNP/Indel: $SNPS_QD/$INDELS_QD, RankSum Min SNP/Indel: $SNPS_RPRS/$INDELS_RPRS, AF Balance: $AF_MIN_BALANCE-$AF_MAX_BALANCE" >> "$ALL_VARIANTS_REPORT"
        echo -e "# CHROM\tPOS\tREF\tALT\tTYPE\tQUAL\tDP\tAD_REF,AD_ALT\tGT\tCALCULATED_AF_ALT\tFINAL_STATUS\tFAIL_REASONS" >> "$ALL_VARIANTS_REPORT"

        # Function to evaluate and report a variant
        report_variant() {
            local TYPE=$1
            local Q_FILTER=$2

            # Filter variants GT=1/1 or GT=0/1 to be relevant, ignoring 0/0.
            # Added || true to the end of the pipe to prevent the script from exiting if there is no output
            bcftools view -H -i 'FMT/GT="1/1" | FMT/GT="0/1"' -v $TYPE "$TEMP_DECOMPOSED" 2>/dev/null | \
            while read -r LINE; do
                # 1. Extract variant data
                CHROM=$(echo "$LINE" | awk '{print $1}')
                POS=$(echo "$LINE" | awk '{print $2}')
                REF=$(echo "$LINE" | awk '{print $4}')
                ALT=$(echo "$LINE" | awk '{print $5}')
                QUAL=$(echo "$LINE" | awk '{print $6}')

                FORMAT_FIELD=$(echo "$LINE" | awk '{print $9}')
                SAMPLE_FIELD=$(echo "$LINE" | awk '{print $10}')

                GT_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^GT$' | cut -d: -f1)
                AD_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^AD$' | cut -d: -f1)
                DP_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^DP$' | cut -d: -f1)

                if [ -z "$GT_INDEX" ] || [ -z "$AD_INDEX" ] || [ -z "$DP_INDEX" ]; then
                    echo -e "$CHROM\t$POS\t$REF\t$ALT\t$TYPE\t$QUAL\tN/A\tN/A\tN/A\tN/A\tPENDING\tMISSING_FORMAT_FIELD" >> "$ALL_VARIANTS_REPORT"
                    continue
                fi

                GT=$(echo "$SAMPLE_FIELD" | cut -d: -f$GT_INDEX)
                AD=$(echo "$SAMPLE_FIELD" | cut -d: -f$AD_INDEX)
                DP=$(echo "$SAMPLE_FIELD" | cut -d: -f$DP_INDEX)
                AD_ALT=$(echo "$AD" | cut -d, -f2)

                AF_CALCULADO="N/A"
                # Added 2>/dev/null to bc to suppress errors if DP is zero/dot
                if [ -n "$DP" ] && [ "$DP" -ne 0 ] && [ "$AD_ALT" != "." ]; then
                    AF_CALCULADO=$(echo "scale=3; $AD_ALT / $DP" | bc 2>/dev/null)
                fi

                # 2. Evaluate FAILURE criteria (to generate REASONS)
                # Added || true to ensure bcftools filter does not cause the script to exit
                IS_Q_PASS=$(bcftools filter -i "${Q_FILTER}" -r $CHROM:$POS-$POS "$TEMP_DECOMPOSED" 2>/dev/null | grep -v '^#' | wc -l || true)

                IS_AF_BALANCE_PASS=1
                if [ "$GT" = "0/1" ]; then
                    # Added || true to ensure bcftools filter does not cause the script to exit
                    IS_AF_BALANCE_PASS=$(bcftools filter -i "${EVAL_AF_BALANCE_FILTER}" -r $CHROM:$POS-$POS "$TEMP_DECOMPOSED" 2>/dev/null | grep -v '^#' | wc -l || true)
                fi

                # 3. Determine FAILURE REASONS
                FINAL_STATUS="PENDING"
                FAIL_REASONS=""

                IS_FAILED_Q=0
                if [ "$IS_Q_PASS" -eq 0 ]; then
                    FAIL_REASONS="${FAIL_REASONS}Q_FAIL(DP_AD_QD_MaxDP);"
                    IS_FAILED_Q=1
                fi

                IS_FAILED_AF=0
                if [ "$GT" = "0/1" ] && [ "$IS_AF_BALANCE_PASS" -eq 0 ]; then
                    FAIL_REASONS="${FAIL_REASONS}AF_IMBALANCE(0.4-0.6);"
                    IS_FAILED_AF=1
                fi

                if [ "$IS_FAILED_Q" -eq 1 ] || [ "$IS_FAILED_AF" -eq 1 ]; then
                    FINAL_STATUS="FAIL_PRELIMINARY"
                fi

                # 4. Write the line to the report
                echo -e "$CHROM\t$POS\t$REF\t$ALT\t$TYPE\t$QUAL\t$DP\t$AD\t$GT\t$AF_CALCULADO\t$FINAL_STATUS\t$FAIL_REASONS" >> "$ALL_VARIANTS_REPORT"

            done || true # CRITICAL: Allows the pipe to fail without exiting the script
        }

        # Generate report for SNPs and then for Indels
        report_variant "snps" "${EVAL_SNPS_Q_FILTER}"
        report_variant "indels" "${EVAL_INDELS_Q_FILTER}"

        REPORT_COUNT_TOTAL=$(wc -l < "$ALL_VARIANTS_REPORT" 2>/dev/null || echo 0)

        if [ "$REPORT_COUNT_TOTAL" -ge 4 ]; then
            REPORT_DATA_LINES=$((REPORT_COUNT_TOTAL - 3))
            echo "    -> ${ALL_VARIANTS_REPORT_SUFFIX} generated with $REPORT_DATA_LINES variants evaluated."
        else
            echo "    -> Warning: ${ALL_VARIANTS_REPORT_SUFFIX} generated, but found few variants to report."
        fi
        # -------------------------------------------------------------------------

        # 2.1.5.b. SNP Filtering (Path A OR Path B)
        echo "        -> 2.1.5.b. Selecting and filtering SNPs..."
        bcftools view -v snps -O z -o "$TEMP_SNPS_VCF" "$TEMP_DECOMPOSED"
        tabix -p vcf "$TEMP_SNPS_VCF"

        bcftools filter -i "${SNPS_FINAL_FILTER}" -o "$TEMP_SNPS_FILTERED" -O z "$TEMP_SNPS_VCF"
        tabix -p vcf "$TEMP_SNPS_FILTERED"

        # 2.1.5.c. Indel Filtering (Path A OR Path B)
        echo "        -> 2.1.5.c. Selecting and filtering Indels..."
        bcftools view -v indels -O z -o "$TEMP_INDELS_VCF" "$TEMP_DECOMPOSED"
        tabix -p vcf "$TEMP_INDELS_VCF"

        bcftools filter -i "${INDELS_FINAL_FILTER}" -o "$TEMP_INDELS_FILTERED" -O z "$TEMP_INDELS_VCF"
        tabix -p vcf "$TEMP_INDELS_FILTERED"

        # 2.1.5.d. Recombination (creates the final PASS VCF)
        echo "        -> 2.1.5.d. Recombining filtered SNPs and Indels. Final output: ${FILTERED_VCF}..."

        VCFS_TO_CONCAT=""
        if [ -s "$TEMP_SNPS_FILTERED" ]; then VCFS_TO_CONCAT="$VCFS_TO_CONCAT $TEMP_SNPS_FILTERED"; fi
        if [ -s "$TEMP_INDELS_FILTERED" ]; then VCFS_TO_CONCAT="$VCFS_TO_CONCAT $TEMP_INDELS_FILTERED"; fi

        if [ -z "$VCFS_TO_CONCAT" ]; then
            echo "        -> Warning: No variants remaining after filtering. Creating valid empty VCF."
            # Creates a valid, but empty, BGZF VCF, keeping the header
            bcftools view -h ${UNPHASED_VCF} 2>/dev/null | bgzip -c > "$FILTERED_VCF"
        else
            bcftools concat $VCFS_TO_CONCAT -a -O z -o "$FILTERED_VCF"
        fi

        # --- 2.1.7. CRITICAL: RECONCILIATION OF TXT REPORT WITH FILTERED VCF ---
        echo "    -> 2.1.7. RECONCILIATION: Correcting PASS/FAIL status in the TXT report..."

        # 1. Extracts the positions (CHROM and POS) of all variants that *actually* passed
        bcftools query -f '%CHROM\t%POS\n' ${FILTERED_VCF} 2>/dev/null > "$PASS_POSITIONS_FILE"

        # 2. Reconciliation with AWK (exact logic maintained for reconciliation)
        TEMP_REPORT="$ALL_VARIANTS_REPORT.tmp"

        awk -v OFS='\t' '
            NR==FNR {
                if (FNR > 0) { pass_positions[$1 ":" $2] = 1; }
                next;
            }
            FNR==1 { print; next; } # Prints the report header

            {
                key = $1 ":" $2;

                if (key in pass_positions) {
                    $11 = "PASS";  # FINAL_STATUS column
                    $12 = "";      # FAIL_REASONS column (cleared)
                }
                else if ($11 == "FAIL_PRELIMINARY") {
                    $11 = "FAIL";
                }
                else if ($11 == "PENDING") {
                    $11 = "FAIL";
                    $12 = "FILTER_FAIL_OTHER";
                }

                print;
            }
        ' "$PASS_POSITIONS_FILE" "$ALL_VARIANTS_REPORT" > "$TEMP_REPORT"

        mv "$TEMP_REPORT" "$ALL_VARIANTS_REPORT"
        echo "    -> Report status corrected: $ALL_VARIANTS_REPORT now matches the final VCF."
        # -------------------------------------------------------------------------


        # 2.1.5.e. CRITICAL: GENERATING UNBALANCED VARIANT REPORT
        echo "        -> 2.1.5.e. Generating Unbalanced Variant Report (0/1 imbalanced)..."

        # COMPLETE EXPRESSION FOR IMBALANCE REPORT
        SNPS_IMBALANCE_REPORT_FILTER="(${NON_RANKSUM_SNPS}) & (${IMBALANCE_AF_CHECK_FILTER})"
        INDELS_IMBALANCE_REPORT_FILTER="(${NON_RANKSUM_INDELS}) & (${IMBALANCE_AF_CHECK_FILTER})"

        # Creates the report header
        echo -e "# Report of High-Quality Unbalanced Variants (FMT/GT=0/1, but AF outside 0.4-0.6)\n# CHROM\tPOS\tREF\t%ALT\tQUAL\tAD\tDP\tGT\tCALCULATED_AF" > "$IMBALANCE_REPORT"

        # Query for unbalanced SNPs (uses the decomposed/normalized VCF)
        # Added || true to the end to prevent the 'set -e' failure if there is no output
        bcftools filter -i "${SNPS_IMBALANCE_REPORT_FILTER}" "$TEMP_DECOMPOSED" 2>/dev/null | \
        bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' >> "$IMBALANCE_REPORT" || true

        # Query for unbalanced INDELs
        # Added || true to the end to prevent the 'set -e' failure if there is no output
        bcftools filter -i "${INDELS_IMBALANCE_REPORT_FILTER}" "$TEMP_DECOMPOSED" 2>/dev/null | \
        bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' >> "$IMBALANCE_REPORT" || true

        REPORT_COUNT=$(wc -l < "$IMBALANCE_REPORT" 2>/dev/null || echo 1)
        echo "        -> ${IMBALANCE_REPORT_SUFFIX} generated with $((REPORT_COUNT - 1)) suspect variant(s)."

    fi

    # --- 2.1.8. DEBUG: List the first 5 variants of each gene/contig ---
    echo "    -> 2.1.8. DEBUG: Listing the first 5 variants detected for EACH contig/gene..."

    # Index the final filtered VCF, as WhatsHap needs it.
    echo "    -> Indexing the filtered VCF for WhatsHap..."
    tabix -p vcf ${FILTERED_VCF}

    # Added || true to find to prevent failure if the VCF is empty (and bcftools finds no contigs)
    ALL_CONTIGS=$(bcftools view -H ${FILTERED_VCF} 2>/dev/null | awk '{print $1}' | sort -u || true)

    if [ -z "$ALL_CONTIGS" ]; then
        echo "        -> Warning: No variants found in the filtered VCF. (Check thresholds)."
    else
        echo "        --- START OF VARIANT REPORT PER CONTIG (MAX 20 LINES EACH) ---"
        for CONTIG in $ALL_CONTIGS; do
            echo ""
            echo "        == CONTIG/GENE: $CONTIG (First 20 Variants) =="
            bcftools query -r $CONTIG -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' ${FILTERED_VCF} | \
            awk 'BEGIN {OFS="\t"} {$6="PASS"; print}' | head -n 20
        done
        echo "        ---------------------------------------------------"
    fi
    # ----------------------------------------------------------------------


    # --- 2.2. PHASING WITH WHATSHAP ---
    echo "    -> 2.2. Phasing variants with WhatsHap (creates $FINAL_VCF_SUFFIX)..."

    whatshap polyphase \
        --ploidy 2 \
        --reference ${GENES_FASTA} \
        --distrust-genotypes \
        ${FILTERED_VCF} \
        ${INPUT_BAM_TO_USE} \
        -o ${PHASED_VCF}

    # Checks if WhatsHap failed
    if [ $? -ne 0 ]; then
        echo "    CRITICAL WHATSHAP ERROR: Phasing failed for isolate $ISOLATE_NAME. Skipping reconstruction."
        cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
        continue
    fi

    # Cleanup of filtered VCF files (which were input for WhatsHap)
    rm -f ${FILTERED_VCF}
    rm -f ${FILTERED_VCF%.*}.tbi

    # --- 2.2.6. INDEXING OF THE PHASED VCF (For the Python script) ---
    echo "    -> 2.2.6. Indexing the phased VCF with tabix..."
    tabix -p vcf ${PHASED_VCF}

    # --- 2.3. HAPLOTYPE RECONSTRUCTION WITH PYTHON ---
    echo "    -> 2.3. Reconstructing haplotypes (Python)..."
    python ${HAPLOTYPE_SCRIPT} --reconstruct_single ${GENES_FASTA} ${RESULTS_DIR} ${ISOLATE_NAME}

    # --- 2.4. FINAL CLEANUP ---
    cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"

    echo "[$CURRENT_COUNT/$NUM_ISOLATES] -> $ISOLATE_NAME COMPLETE."
    echo "" # Empty line for clarity
done

# ---- PHASE 3: Final Consolidation ----

echo "### STEP 3: Final Consolidation"
echo "------------------------------------------------------"
echo "3. Consolidating all haplotypes into a FASTA per gene..."
echo "------------------------------------------------------"

# Added || true in case the python script fails due to lack of output (but output should be generated)
python ${HAPLOTYPE_SCRIPT} --consolidate_all ${RESULTS_DIR} || true

echo "------------------------------------------------------"
echo "PIPELINE COMPLETE! Final FASTA files per gene are in $RESULTS_DIR."
echo "------------------------------------------------------"
