# HapHap: Automated Phasing and Protein-Level Variant Discovery Pipeline

## Overview

**HapHap** is a comprehensive pipeline for analyzing diploid sequencing data and reconstructing phased haplotypes per gene. It combines robust read mapping, variant calling, filtering, phasing, and final protein-level mutation annotation. The workflow is split into two main components:

1. **Phasing and Haplotype Pipeline (`phasing_pipeline.sh`):**
    - Automates mapping (`bwa` + `samtools`), variant calling (`GATK`), and intricate filtering (`bcftools`), followed by read-backed phasing with `WhatsHap`.
    - Produces phased, quality-filtered VCFs and generates FASTA files of reconstructed sequences per haplotype for each gene/isolate.

2. **Haplotype Mutation Annotation (`haplotype_reconstruction.py`):**
    - Given the pipeline output (per-gene haplotype FASTAs), and a reference CDS FASTA, identifies all amino acid (AA) substitutions, in-frame indels, and frameshift events, reporting both zygosity (hom/het) and primary sequence alignments.
    - Generates ready-to-analyze CSV tables, detailed protein alignment logs, and summary reports for downstream analysis.

These programs are fully interoperable and designed for fungal population genomics and clinical isolate comparative studies.

---

## Features

- **Automated, error-resilient pipeline:** Handles multiple samples with strict reporting and intermediate quality checks.
- **Gene-level and sample-level flexibility:** Restricts analysis to specified genes or all loci in the reference FASTA.
- **Dynamic, data-driven filtering:** Depth and quality thresholds calculated per sample for robust variant retention.
- **Comprehensive mutation profiling:** Final Python analysis calls over all phased haplotypes with explicit frameshift detection and annotation.
- **Modularity:** Core scripts can be run independently or as a two-step workflow (see usage below).

---

## Dependencies

### `phasing_pipeline.sh` (Step 1)

Requires:

- **BWA** (read mapping)
- **SAMtools** (BAM manipulation)
- **GATK** (variant calling, reference indexing)
- **bcftools** (variant filtering, normalization)
- **WhatsHap** (read-backed genotype phasing)
- **Python3**
- **haplotype_reconstruction.py** (called internally for reconstruction steps)

### `haplotype_reconstruction.py` (Step 2)

Requires Python libraries:

- `biopython` (>=1.76)
- `pandas`
- `tqdm`
- (No external compiled dependencies—uses pure Python and system libraries.)

Install via:
```sh
pip install biopython pandas tqdm
```

---

## Installation & Setup

1. **Clone the repository:**

```sh
git clone https://github.com/zagalom/HapHap.git
cd HapHap
```

2. **Install dependencies** (see above).

3. **Prepare reference files and sample FASTQ data:**
   - Concatenated reference genes FASTA (`genes.fasta`)
   - Folder of paired-end FASTQs, named as `SAMPLE_1.fastq.gz`, `SAMPLE_2.fastq.gz`
   - Optionally, a comma-separated list of genes (e.g., `ERG11,GSC1,UPC2`)

---

## Usage

### 1. Run the Phasing Pipeline

```sh
bash phasing_pipeline.sh <genes_fasta> <results_dir> <fastq_dir> [threads] [genes_list]
```

- `genes_fasta` : Reference genes FASTA file (required)
- `results_dir` : Output directory (will contain per-sample subfolders)
- `fastq_dir`   : Directory containing paired FASTQs (required)
- `[threads]`   : Number of threads (default: 4, optional)
- `[genes_list]`: Comma-separated list of gene names (optional, restricts analysis)

**Example:**
```sh
bash phasing_pipeline.sh ref_genes.fasta HapHap_results samples/ 8 ERG11,GSC1
```

This will create phased VCFs and FASTA files per sample and gene in the output directory.

### 2. Run Protein-Variant Annotation

After the first step, run the Python script to annotate all mutations by gene:

```sh
python haplotype_reconstruction.py <haplotype_fasta_dir> <reference_genes_fasta> <output_dir>
```

- `haplotype_fasta_dir` : Directory of phased haplotype FASTA files (output from previous step)
- `reference_genes_fasta`: Reference CDS FASTA file (same as used above)
- `output_dir`: Directory to save CSV tables and reports

**Example:**
```sh
python haplotype_reconstruction.py HapHap_results/ ref_genes.fasta HapHap_annotation_output/
```

---

## Workflow Integration

The recommended workflow is:

1. **Run `phasing_pipeline.sh`** to produce sample/gene-level phased haplotype FASTAs.
2. **Run `haplotype_reconstruction.py`** to annotate every variant (AA, indel, frameshift) in each isolate/gene against the reference.

The shell script automatically calls the Python module for basic FASTA consolidation, but for full mutation profiling, run the Python script as described above.

---

## Expected Output

- **Phased & Filtered VCFs**: Per isolate, per gene.
- **Haplotype FASTAs**: Per gene, including both reconstructed haplotypes per sample.
- **Variant CSV Tables**: For each gene, shows all isolates and their detected mutations, marked as `Hom` or `Het`.
- **Summary Reports & Alignment Logs**: For each gene, detailed protein alignments and variant listings.

---

## Notes & Best Practices

- Always verify that your sample FASTQ names follow expected naming (`SAMPLE_1.fastq.gz`, `SAMPLE_2.fastq.gz`).
- The pipeline is robust to missing pairs, low-depth, or suboptimal samples—it will skip and report problematic isolates.
- For custom filtering thresholds (e.g., adjusting DP/AD/QD), manually edit the shell script variables.

---

## Citation & Contact

Developed by [@zagalom](https://github.com/zagalom)  
For questions, bug reports, or collaborations, please open [an issue](https://github.com/zagalom/HapHap/issues).

---

## License

Released under the MIT license. See [LICENSE](LICENSE) for details.
