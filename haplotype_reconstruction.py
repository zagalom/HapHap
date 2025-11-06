#!/usr/bin/env python3
import sys
import os
import re
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------

def load_reference_genes(reference_fasta):
    """Return a list of contig IDs (gene names) from the reference FASTA."""
    genes = []
    try:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            genes.append(record.id)
    except FileNotFoundError:
        sys.exit(f"Fatal Error: Reference FASTA file not found at {reference_fasta}")
    return genes

def parse_gene_list(reference_fasta, genes_arg):
    """Return list of target genes: from user or from FASTA if not provided."""
    if genes_arg:
        genes = [g.strip() for g in genes_arg.split(",") if g.strip()]
        print(f"    -> Genes manually specified: {', '.join(genes)}")
    else:
        # Load all genes from the reference FASTA
        genes = load_reference_genes(reference_fasta)
        print(f"    -> No gene specified: all {len(genes)} genes from FASTA will be used.")
    return genes

# -------------------------------------------------------------
# Main modes
# -------------------------------------------------------------

def reconstruct_single(reference_fasta, results_dir, isolate_name, genes):
    print(f"\n🚀 Reconstructing haplotypes for isolate {isolate_name}...")
    vcf_path = os.path.join(results_dir, isolate_name, f"{isolate_name}.filtered.phased.vcf.gz")

    if not os.path.exists(vcf_path):
        print(f"    ⚠️ VCF file not found: {vcf_path}")
        return

    fasta_records = []
    vcf_in = pysam.VariantFile(vcf_path)
    
    # Load reference sequences into a dictionary for quick access
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
    except FileNotFoundError:
        print(f"    ⚠️ Reference FASTA file not found: {reference_fasta}")
        return

    for gene in genes:
        if gene not in vcf_in.header.contigs:
            print(f"    ⚠️ Gene {gene} not found in VCF header, skipping.")
            continue
            
        if gene not in ref_dict:
            print(f"    ⚠️ Gene {gene} not found in reference FASTA.")
            continue

        gene_seq = ref_dict[gene].seq

        # Start with the reference sequence as two lists of characters
        hap1, hap2 = list(gene_seq), list(gene_seq)
        
        # Keep track of insertion/deletion offsets for each haplotype
        offsets = {0: 0, 1: 0} # offset for hap1 (GT[0]), offset for hap2 (GT[1])

        for rec in vcf_in.fetch(gene):
            if rec.pos <= 0:
                continue
            
            # Ensure we have a phased genotype (GT)
            if "GT" not in rec.samples[0] or \
               rec.samples[0]["GT"] is None or \
               len(rec.samples[0]["GT"]) < 2 or \
               rec.samples[0]["GT"][0] is None or \
               rec.samples[0]["GT"][1] is None:
                continue
                
            gt = rec.samples[0]["GT"]
            alleles = rec.alleles
            ref_allele = alleles[0]
            
            # Process Haplotype 1 (gt[0]) and Haplotype 2 (gt[1])
            for i, hap_index in enumerate([gt[0], gt[1]]):
                hap_list = hap1 if i == 0 else hap2
                offset = offsets[i]
                
                alt_allele = alleles[hap_index]
                
                start_pos = rec.pos - 1 # 0-based
                end_pos = start_pos + len(ref_allele)
                
                # Apply the variant based on allele lengths
                if len(ref_allele) == len(alt_allele): # SNP or MNV
                    # Replace the segment
                    hap_list[start_pos + offset : end_pos + offset] = list(alt_allele)
                
                elif len(ref_allele) > len(alt_allele): # Deletion (or complex indel starting with deletion)
                    # Remove the reference segment
                    del hap_list[start_pos + offset : end_pos + offset]
                    # Insert the alternative sequence (can be an empty string or a short sequence)
                    hap_list[start_pos + offset : start_pos + offset] = list(alt_allele)
                    # Update the offset (the sequence got shorter)
                    offsets[i] -= (len(ref_allele) - len(alt_allele))
                    
                else: # Insertion (or complex indel starting with insertion)
                    # Remove the reference sequence (usually 1 base)
                    del hap_list[start_pos + offset : end_pos + offset]
                    # Insert the alternative sequence
                    hap_list[start_pos + offset : start_pos + offset] = list(alt_allele)
                    # Update the offset (the sequence got longer)
                    offsets[i] += (len(alt_allele) - len(ref_allele))

        # Create SeqRecord objects
        fasta_records.append(SeqRecord(Seq("".join(hap1)), id=f"{gene}_Haplotype1_{isolate_name}", description=""))
        fasta_records.append(SeqRecord(Seq("".join(hap2)), id=f"{gene}_Haplotype2_{isolate_name}", description=""))

    output_fasta = os.path.join(results_dir, isolate_name, f"{isolate_name}_haplotypes.fasta")
    SeqIO.write(fasta_records, output_fasta, "fasta")
    print(f"    ✅ File generated: {output_fasta}")

# -------------------------------------------------------------
# Consolidate mode
# -------------------------------------------------------------

def consolidate_all(results_dir, genes):
    print("\n📦 Consolidating haplotypes from all isolates...")
    fasta_files = []
    # Recursively find all *_haplotypes.fasta files in subdirectories
    for root, _, files in os.walk(results_dir):
        for f in files:
            if f.endswith("_haplotypes.fasta"):
                fasta_files.append(os.path.join(root, f))

    if not fasta_files:
        print("⚠️ No haplotype files (*_haplotypes.fasta) found.")
        return

    print(f"    -> Found {len(fasta_files)} _haplotypes.fasta files to process.")

    for gene in genes:
        records = []
        for fpath in fasta_files:
            try:
                # Iterate through records in the isolate's haplotype file
                for rec in SeqIO.parse(fpath, "fasta"):
                    # Only append records matching the current gene
                    if rec.id.startswith(gene + "_"):
                        records.append(rec)
            except Exception as e:
                print(f"    ⚠️ Error reading {fpath}: {e}")
                
        if not records:
            print(f"    ⚠️ No haplotype found for {gene} in the processed files.")
            continue
            
        # Output filename MUST be GENE.fasta for the subsequent analysis script
        output_path = os.path.join(results_dir, f"{gene}.fasta")
        SeqIO.write(records, output_path, "fasta")
        print(f"    ✅ {gene}: {len(records)} haplotypes written -> {output_path}")

# -------------------------------------------------------------
# Entry point
# -------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python haplotype_reconstruction.py [--reconstruct_single|--consolidate_all] <args> [--genes gene1,gene2,...]")

    # Parse optional --genes argument
    genes_arg = None
    if "--genes" in sys.argv:
        try:
            idx = sys.argv.index("--genes")
            # Extract the gene list and remove --genes and its argument from argv
            genes_arg = sys.argv[idx + 1]
            sys.argv = sys.argv[:idx] + sys.argv[idx+2:] 
        except IndexError:
            sys.exit("Error: --genes requires an argument (e.g., ERG11,GSC1)")
        except ValueError:
            pass # --genes not found

    if sys.argv[1] == "--reconstruct_single" and len(sys.argv) >= 5:
        # Expected args: reference_fasta, results_dir, isolate_name
        ref, resdir, isolate = sys.argv[2:5]
        genes = parse_gene_list(ref, genes_arg)
        reconstruct_single(ref, resdir, isolate, genes)
        
    elif sys.argv[1] == "--consolidate_all" and len(sys.argv) >= 4:
        # Expected args: results_dir, reference_fasta
        resdir = sys.argv[2]
        fasta_path = sys.argv[3] # Path to the reference FASTA
        
        if not os.path.exists(fasta_path):
            sys.exit(f"Error: Reference FASTA file not found at {fasta_path}")

        genes = parse_gene_list(fasta_path, genes_arg)
        consolidate_all(resdir, genes)
        
    else:
        print("Error: Invalid arguments for the selected mode.")
        print("Mode --reconstruct_single requires: <reference_fasta> <results_dir> <isolate_name>")
        print("Mode --consolidate_all requires: <results_dir> <reference_fasta>")
        sys.exit(1)
