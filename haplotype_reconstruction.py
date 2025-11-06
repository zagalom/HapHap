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
    """Return a list of contig IDs from the reference FASTA."""
    genes = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        genes.append(record.id)
    return genes

def parse_gene_list(reference_fasta, genes_arg):
    """Return list of target genes: from user or from FASTA if not provided."""
    if genes_arg:
        genes = [g.strip() for g in genes_arg.split(",") if g.strip()]
        print(f"   -> Genes especificados manualmente: {', '.join(genes)}")
    else:
        genes = load_reference_genes(reference_fasta)
        print(f"   -> Nenhum gene especificado: serão usados todos os {len(genes)} genes do FASTA.")
    return genes

# -------------------------------------------------------------
# Main modes
# -------------------------------------------------------------

def reconstruct_single(reference_fasta, results_dir, isolate_name, genes):
    print(f"\n🚀 A reconstruir haplótipos para o isolado {isolate_name}...")
    vcf_path = os.path.join(results_dir, isolate_name, f"{isolate_name}.filtered.phased.vcf.gz")

    if not os.path.exists(vcf_path):
        print(f"   ⚠️ Ficheiro VCF não encontrado: {vcf_path}")
        return

    fasta_records = []
    vcf_in = pysam.VariantFile(vcf_path)

    for gene in genes:
        if gene not in vcf_in.header.contigs:
            print(f"   ⚠️ Gene {gene} não encontrado no VCF, ignorado.")
            continue
        try:
            ref_seq = next(SeqIO.parse(reference_fasta, "fasta"))
        except Exception:
            ref_seq = None
        gene_seq = next((r.seq for r in SeqIO.parse(reference_fasta, "fasta") if r.id == gene), None)
        if not gene_seq:
            print(f"   ⚠️ Gene {gene} não encontrado no FASTA.")
            continue

        hap1, hap2 = list(gene_seq), list(gene_seq)
        offset1 = offset2 = 0

        for rec in vcf_in.fetch(gene):
            if rec.pos <= 0:
                continue
            alleles = rec.alleles
            gt = rec.samples[0]["GT"]
            if not gt or len(gt) < 2:
                continue

            ref_allele = alleles[0]
            alt1 = alleles[gt[0]] if gt[0] is not None else ref_allele
            alt2 = alleles[gt[1]] if gt[1] is not None else ref_allele

            start = rec.pos - 1
            end = start + len(ref_allele)

            for hap, alt, offset in [(hap1, alt1, offset1), (hap2, alt2, offset2)]:
                if len(ref_allele) == len(alt):  # SNP
                    hap[start + offset:end + offset] = list(alt)
                elif len(ref_allele) > len(alt):  # deletion
                    del hap[start + offset:end + offset]
                    hap[start + offset:start + offset] = list(alt)
                    offset -= len(ref_allele) - len(alt)
                else:  # insertion
                    hap[start + offset:start + offset + len(ref_allele)] = list(alt)
                    offset += len(alt) - len(ref_allele)

        fasta_records.append(SeqRecord(Seq("".join(hap1)), id=f"{gene}_Haplotype1_{isolate_name}", description=""))
        fasta_records.append(SeqRecord(Seq("".join(hap2)), id=f"{gene}_Haplotype2_{isolate_name}", description=""))

    output_fasta = os.path.join(results_dir, isolate_name, f"{isolate_name}_haplotypes.fasta")
    SeqIO.write(fasta_records, output_fasta, "fasta")
    print(f"   ✅ Ficheiro gerado: {output_fasta}")

# -------------------------------------------------------------
# Consolidate mode
# -------------------------------------------------------------

def consolidate_all(results_dir, genes):
    print("\n📦 A consolidar haplótipos de todos os isolados...")
    fasta_files = []
    for root, _, files in os.walk(results_dir):
        for f in files:
            if f.endswith("_haplotypes.fasta"):
                fasta_files.append(os.path.join(root, f))

    if not fasta_files:
        print("⚠️ Nenhum ficheiro de haplótipos encontrado.")
        return

    for gene in genes:
        records = []
        for fpath in fasta_files:
            for rec in SeqIO.parse(fpath, "fasta"):
                if rec.id.startswith(gene + "_"):
                    records.append(rec)
        if not records:
            print(f"   ⚠️ Nenhum haplótipo encontrado para {gene}")
            continue
        output_path = os.path.join(results_dir, f"{gene}_all_isolates.fasta")
        SeqIO.write(records, output_path, "fasta")
        print(f"   ✅ {gene}: {len(records)} haplótipos gravados -> {output_path}")

# -------------------------------------------------------------
# Entry point
# -------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Uso: python haplotype_reconstruction.py [--reconstruct_single|--consolidate_all] <args> [--genes gene1,gene2,...]")

    # Parse optional --genes argument
    genes_arg = None
    if "--genes" in sys.argv:
        idx = sys.argv.index("--genes")
        genes_arg = sys.argv[idx + 1]
        sys.argv = sys.argv[:idx]  # remove from argument list

    if sys.argv[1] == "--reconstruct_single" and len(sys.argv) >= 5:
        ref, resdir, isolate = sys.argv[2:5]
        genes = parse_gene_list(ref, genes_arg)
        reconstruct_single(ref, resdir, isolate, genes)
    elif sys.argv[1] == "--consolidate_all" and len(sys.argv) >= 3:
        resdir = sys.argv[2]
        # fallback reference fasta for listing contigs (optional if user gave genes)
        fasta_path = os.path.join(resdir, "reference.fasta") if os.path.exists(os.path.join(resdir, "reference.fasta")) else None
        genes = parse_gene_list(fasta_path, genes_arg) if fasta_path else [g for g in (genes_arg or "").split(",") if g]
        consolidate_all(resdir, genes)
    else:
        sys.exit("Erro: argumentos inválidos.")
