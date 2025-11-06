#!/bin/bash

set -euo pipefail

# ---- CONFIGURATION ----

# Allow user to provide input or prompt if not supplied
GENES_FASTA=${1:-}
RESULTS_DIR=${2:-}
FASTQ_DIR=${3:-}
THREADS=${4:-4}  # Default to 4 if threads not specified

function print_usage {
    echo "Usage:"
    echo "  $0 <genes_fasta> <results_dir> <fastq_dir> [threads]"
    echo "    <genes_fasta>  : Reference genes FASTA"
    echo "    <results_dir>  : Directory for output/sample subdirectories"
    echo "    <fastq_dir>    : Directory containing paired FASTQs (SAMPLE_1.fastq.gz/SAMPLE_2.fastq.gz)"
    echo "    [threads]      : Number of threads to use (default: 4)"
    exit 1
}

[ -n "$GENES_FASTA" ] && [ -n "$RESULTS_DIR" ] && [ -n "$FASTQ_DIR" ] || print_usage

HAPLOTYPE_SCRIPT="haplotype_reconstruction.py"
export LC_NUMERIC="C"   # Ensure decimal points for awk/bc

DP_MULTIPLIER="2"
BYPASS_FILTER="false"

# For GATK Read Group info (update if necessary)
RG_ID="A"
RG_LB="A"
RG_PL="ILLUMINA"
RG_PU="L001"

TEMP_VCF_SUFFIX=".unphased.vcf.gz"
FILTERED_VCF_SUFFIX=".filtered.vcf.gz"
FINAL_VCF_SUFFIX=".phased.vcf.gz"
IMBALANCE_REPORT_SUFFIX=".imbalanced_report.txt"
ALL_VARIANTS_REPORT_SUFFIX=".all_variants_filter_report.txt"

# ---- PHASE 0: Mapping FASTQs and generating BAMs ----

echo "### STEP 0: Mapping reads and generating sorted BAMs..."

if ! ls "${GENES_FASTA}".*bwt &>/dev/null; then
    echo "Indexing $GENES_FASTA for BWA..."
    bwa index "$GENES_FASTA"
fi

mkdir -p "$RESULTS_DIR"

for fq1 in "$FASTQ_DIR"/*_1*.fastq*; do
    [ -e "$fq1" ] || continue
    SAMPLE=$(basename "$fq1" | sed 's/_1.*//')
    fq2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    [ -f "$fq2" ] || fq2="$FASTQ_DIR/${SAMPLE}_2.fastq"
    if [ ! -f "$fq2" ]; then
        echo "Warning: Missing 2 file for $fq1. Skipping sample $SAMPLE."
        continue
    fi
    smdir="$RESULTS_DIR/$SAMPLE"
    mkdir -p "$smdir"
    outbam="$smdir/${SAMPLE}.sorted.bam"
    if [ -s "$outbam" ]; then
        echo "[$SAMPLE] sorted BAM exists. Skipping alignment."
        continue
    fi
    echo "[$SAMPLE] Aligning and sorting using BWA/SAMtools..."
    bwa mem -t "$THREADS" "$GENES_FASTA" "$fq1" "$fq2" |
        samtools view -bh -@ "$THREADS" - |
        samtools sort -@ "$THREADS" -o "$outbam"
    samtools index "$outbam"
    echo "[$SAMPLE] Finished BAM: $outbam"
done

# ---- Prepare reference indices for variant calling ----

echo "### STEP 1: Preparing reference indices..."

if [ ! -f "${GENES_FASTA%.*}.dict" ]; then
    echo "Creating GATK sequence dictionary for reference..."
    gatk CreateSequenceDictionary -R "$GENES_FASTA"
fi

if [ ! -f "${GENES_FASTA}.fai" ]; then
    echo "Creating FASTA index (.fai) for $GENES_FASTA..."
    samtools faidx "$GENES_FASTA"
fi

# ---- PHASE 1: Per-sample variant calling, phasing, and haplotype reconstruction ----

echo "### STEP 2: Sample-wise variant calling, phasing, haplotype reconstruction..."

ALL_SAMPLES=$(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
NUM_SAMPLES=$(echo "$ALL_SAMPLES" | wc -l)
if [ "$NUM_SAMPLES" -eq 0 ]; then
    echo "ERROR: No sample folders found in $RESULTS_DIR"
    exit 1
fi

CURRENT=0
for SAMPLE_DIR in $ALL_SAMPLES; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    ((CURRENT++))
    echo
    echo "------ [$CURRENT/$NUM_SAMPLES] SAMPLE: $SAMPLE ------"
    BAM_FILE="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
    RG_BAM="$SAMPLE_DIR/${SAMPLE}.rg.bam"
    UNPHASED_VCF="$SAMPLE_DIR/${SAMPLE}${TEMP_VCF_SUFFIX}"
    FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}${FILTERED_VCF_SUFFIX}"
    PHASED_VCF="$SAMPLE_DIR/${SAMPLE}${FINAL_VCF_SUFFIX}"
    IMBALANCE_REPORT="$SAMPLE_DIR/${SAMPLE}${IMBALANCE_REPORT_SUFFIX}"
    ALL_VARIANTS_REPORT="$SAMPLE_DIR/${SAMPLE}${ALL_VARIANTS_REPORT_SUFFIX}"

    # Clean up temp files for reruns
    rm -f "$RG_BAM" "$RG_BAM.bai" "$UNPHASED_VCF" "$UNPHASED_VCF.tbi"

    # 1. Add/replace read groups
    if [ ! -f "$BAM_FILE" ]; then
        echo "Warning: BAM not found for $SAMPLE, skipping."
        continue
    fi
    echo "  [1] Updating BAM read groups..."
    gatk AddOrReplaceReadGroups \
        -I "$BAM_FILE" -O "$RG_BAM" \
        -RGID "$RG_ID" -RGLB "$RG_LB" -RGPL "$RG_PL" -RGPU "$RG_PU" -RGSM "$SAMPLE" --CREATE_INDEX true

    # 2. Variant calling
    echo "  [2] Variant calling (GATK HaplotypeCaller)..."
    gatk HaplotypeCaller -R "$GENES_FASTA" -I "$RG_BAM" -O "$UNPHASED_VCF" \
        --minimum-mapping-quality 30 -ploidy 2 --sample-name "$SAMPLE" --min-base-quality-score 20

    # 3. Dynamic depth filter calculation
    echo "  [3] Calculating dynamic DP max threshold..."
    DP_VALUES=$(bcftools query -f '[%DP]\n' "$UNPHASED_VCF" -i 'FORMAT/DP>0' 2>/dev/null)
    NUM_VARIANTS=$(echo "$DP_VALUES" | wc -l)
    SNPS_DP_TOTAL=20
    MAX_DP_FILTER=1000
    if [ "$NUM_VARIANTS" -gt 0 ]; then
        DP_MEAN=$(echo "$DP_VALUES" | awk 'BEGIN {sum=0} {sum+=$1} END {printf "%.3f", sum/NR}')
        NEW_MAX_DP=$(echo "$DP_MEAN * $DP_MULTIPLIER" | bc)
        NEW_MAX_DP_ROUND=$(printf "%.0f" "$NEW_MAX_DP")
        MIN_LIMIT=$((SNPS_DP_TOTAL * 2))
        if [ "$NEW_MAX_DP_ROUND" -ge "$MIN_LIMIT" ]; then
            MAX_DP_FILTER="$NEW_MAX_DP_ROUND"
        fi
        echo "    Mean DP: $DP_MEAN, Max DP filter: $MAX_DP_FILTER"
    fi

    # 4. Filtering
    SNPS_RPRS="-8.0"
    SNPS_BQRS="-8.0"
    SNPS_AD_ALT="10"
    SNPS_QD="2.0"
    INDELS_RPRS="-20.0"
    INDELS_BQRS="-8.0"
    INDELS_AD_ALT="10"
    INDELS_QD="10.0"
    AF_MIN_BALANCE="0.4"
    AF_MAX_BALANCE="0.6"
    NON_RANKSUM_SNPS="FMT/AD[0:1] > ${SNPS_AD_ALT} & FMT/DP[0] > ${SNPS_DP_TOTAL} & FMT/DP[0] < ${MAX_DP_FILTER} & INFO/QD > ${SNPS_QD}"
    NON_RANKSUM_INDELS="FMT/AD[0:1] > ${INDELS_AD_ALT} & FMT/DP[0] > $SNPS_DP_TOTAL & FMT/DP[0] < $MAX_DP_FILTER & INFO/QD > ${INDELS_QD}"
    PERFECT_1_1_SNPS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & ($NON_RANKSUM_SNPS)"
    PERFECT_1_1_INDELS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & ($NON_RANKSUM_INDELS)"
    STRICT_Q_SNPS="(ReadPosRankSum > $SNPS_RPRS & BaseQRankSum > $SNPS_BQRS)"
    STRICT_Q_INDELS="(ReadPosRankSum > $INDELS_RPRS & BaseQRankSum > $INDELS_BQRS)"
    STANDARD_AF_GT_CHECK="(FMT/GT=\"1/1\" | (FMT/GT=\"0/1\" & (FMT/AD[0:1] / FMT/DP[0]) >= $AF_MIN_BALANCE & (FMT/AD[0:1] / FMT/DP[0]) <= $AF_MAX_BALANCE))"
    STANDARD_PASS_SNPS="($NON_RANKSUM_SNPS & $STRICT_Q_SNPS & $STANDARD_AF_GT_CHECK)"
    STANDARD_PASS_INDELS="($NON_RANKSUM_INDELS & $STRICT_Q_INDELS & $STANDARD_AF_GT_CHECK)"
    SNPS_FINAL_FILTER="($PERFECT_1_1_SNPS) | ($STANDARD_PASS_SNPS)"
    INDELS_FINAL_FILTER="($PERFECT_1_1_INDELS) | ($STANDARD_PASS_INDELS)"

    echo "  [4] Filtering and phasing..."
    bcftools norm -f "$GENES_FASTA" -m -any -O v -o "$SAMPLE_DIR/decomposed.vcf" "$UNPHASED_VCF"
    bcftools view -v snps -O z -o "$SAMPLE_DIR/snps.vcf.gz" "$SAMPLE_DIR/decomposed.vcf"
    bcftools view -v indels -O z -o "$SAMPLE_DIR/indels.vcf.gz" "$SAMPLE_DIR/decomposed.vcf"
    tabix -p vcf "$SAMPLE_DIR/snps.vcf.gz"
    tabix -p vcf "$SAMPLE_DIR/indels.vcf.gz"
    bcftools filter -i "$SNPS_FINAL_FILTER" -O z -o "$SAMPLE_DIR/snps.filtered.vcf.gz" "$SAMPLE_DIR/snps.vcf.gz"
    bcftools filter -i "$INDELS_FINAL_FILTER" -O z -o "$SAMPLE_DIR/indels.filtered.vcf.gz" "$SAMPLE_DIR/indels.vcf.gz"
    tabix -p vcf "$SAMPLE_DIR/snps.filtered.vcf.gz"
    tabix -p vcf "$SAMPLE_DIR/indels.filtered.vcf.gz"
    bcftools concat "$SAMPLE_DIR/snps.filtered.vcf.gz" "$SAMPLE_DIR/indels.filtered.vcf.gz" -a -O z -o "$FILTERED_VCF" || \
        (bcftools view -h "$UNPHASED_VCF" | bgzip -c > "$FILTERED_VCF")
    tabix -p vcf "$FILTERED_VCF"

    # 5. Phasing with WhatsHap
    echo "  [5] Variant phasing with WhatsHap..."
    whatshap polyphase --ploidy 2 --reference "$GENES_FASTA" --distrust-genotypes "$FILTERED_VCF" "$RG_BAM" -o "$PHASED_VCF"

    tabix -p vcf "$PHASED_VCF"

    # 6. Haplotype reconstruction
    echo "  [6] Reconstructing haplotypes with python..."
    python "$HAPLOTYPE_SCRIPT" --reconstruct_single "$GENES_FASTA" "$RESULTS_DIR" "$SAMPLE"

    # Clean up sample temp files to save space
    rm -f "$RG_BAM" "$RG_BAM.bai" "$UNPHASED_VCF" "$UNPHASED_VCF.tbi" \
          "$SAMPLE_DIR/"*.vcf "$SAMPLE_DIR/"*.vcf.gz "$SAMPLE_DIR/"*.tbi 2>/dev/null || true
    echo "[FINISHED] $SAMPLE"
done

# ---- FINAL CONSOLIDATION ----

echo
echo "### STEP 3: Consolidating all haplotypes per gene..."
python "$HAPLOTYPE_SCRIPT" --consolidate_all "$RESULTS_DIR"

echo
echo "### PIPELINE COMPLETE ###"
echo "Output is in $RESULTS_DIR. Please check per-sample '.imbalanced_report.txt' and '.all_variants_filter_report.txt' for possible variant calls quality issues."
