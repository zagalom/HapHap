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
    try:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            genes.append(record.id)
    except FileNotFoundError:
        sys.exit(f"Erro fatal: Ficheiro de referência FASTA não encontrado em {reference_fasta}")
    return genes

def parse_gene_list(reference_fasta, genes_arg):
    """Return list of target genes: from user or from FASTA if not provided."""
    if genes_arg:
        genes = [g.strip() for g in genes_arg.split(",") if g.strip()]
        print(f"    -> Genes especificados manualmente: {', '.join(genes)}")
    else:
        genes = load_reference_genes(reference_fasta)
        print(f"    -> Nenhum gene especificado: serão usados todos os {len(genes)} genes do FASTA.")
    return genes

# -------------------------------------------------------------
# Main modes
# -------------------------------------------------------------

def reconstruct_single(reference_fasta, results_dir, isolate_name, genes):
    print(f"\n🚀 A reconstruir haplótipos para o isolado {isolate_name}...")
    vcf_path = os.path.join(results_dir, isolate_name, f"{isolate_name}.filtered.phased.vcf.gz")

    if not os.path.exists(vcf_path):
        print(f"    ⚠️ Ficheiro VCF não encontrado: {vcf_path}")
        return

    fasta_records = []
    vcf_in = pysam.VariantFile(vcf_path)

    # Carregar referências num dicionário para acesso rápido
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
    except FileNotFoundError:
        print(f"    ⚠️ Ficheiro de referência FASTA não encontrado: {reference_fasta}")
        return

    for gene in genes:
        if gene not in vcf_in.header.contigs:
            print(f"    ⚠️ Gene {gene} não encontrado no cabeçalho VCF, ignorado.")
            continue

        if gene not in ref_dict:
            print(f"    ⚠️ Gene {gene} não encontrado no FASTA de referência.")
            continue

        gene_seq = ref_dict[gene].seq

        hap1, hap2 = list(gene_seq), list(gene_seq)
        offset1, offset2 = 0, 0

        # Variáveis para armazenar os offsets
        offsets = {0: 0, 1: 0} # offset para hap1 (GT[0]), offset para hap2 (GT[1])

        for rec in vcf_in.fetch(gene):
            if rec.pos <= 0:
                continue

            # Garantir que temos um genótipo fasedo
            if "GT" not in rec.samples[0] or \
               rec.samples[0]["GT"] is None or \
               len(rec.samples[0]["GT"]) < 2 or \
               rec.samples[0]["GT"][0] is None or \
               rec.samples[0]["GT"][1] is None:
                continue

            gt = rec.samples[0]["GT"]
            alleles = rec.alleles
            ref_allele = alleles[0]

            # Processar Haplotype 1 (gt[0]) e Haplotype 2 (gt[1])
            for i, hap_index in enumerate([gt[0], gt[1]]):
                hap_list = hap1 if i == 0 else hap2
                offset = offsets[i]

                alt_allele = alleles[hap_index]

                start_pos = rec.pos - 1 # 0-based
                end_pos = start_pos + len(ref_allele)

                # Aplicar a variante
                if len(ref_allele) == len(alt_allele): # SNP ou MNV
                    hap_list[start_pos + offset : end_pos + offset] = list(alt_allele)

                elif len(ref_allele) > len(alt_allele): # Deletion
                    # Remover a sequência de referência
                    del hap_list[start_pos + offset : end_pos + offset]
                    # Inserir a sequência alternativa (que pode ser 1 base, ex: 'TA' -> 'T')
                    hap_list[start_pos + offset : start_pos + offset] = list(alt_allele)
                    # Atualizar o offset
                    offsets[i] -= (len(ref_allele) - len(alt_allele))

                else: # Insertion
                    # Remover a sequência de referência (normalmente 1 base)
                    del hap_list[start_pos + offset : end_pos + offset]
                    # Inserir a sequência alternativa
                    hap_list[start_pos + offset : start_pos + offset] = list(alt_allele)
                    # Atualizar o offset
                    offsets[i] += (len(alt_allele) - len(ref_allele))

        fasta_records.append(SeqRecord(Seq("".join(hap1)), id=f"{gene}_Haplotype1_{isolate_name}", description=""))
        fasta_records.append(SeqRecord(Seq("".join(hap2)), id=f"{gene}_Haplotype2_{isolate_name}", description=""))

    output_fasta = os.path.join(results_dir, isolate_name, f"{isolate_name}_haplotypes.fasta")
    SeqIO.write(fasta_records, output_fasta, "fasta")
    print(f"    ✅ Ficheiro gerado: {output_fasta}")

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
        print("⚠️ Nenhum ficheiro de haplótipos (*_haplotypes.fasta) foi encontrado.")
        return

    print(f"    -> Encontrados {len(fasta_files)} ficheiros _haplotypes.fasta para processar.")

    for gene in genes:
        records = []
        for fpath in fasta_files:
            try:
                for rec in SeqIO.parse(fpath, "fasta"):
                    if rec.id.startswith(gene + "_"):
                        records.append(rec)
            except Exception as e:
                print(f"    ⚠️ Erro ao ler {fpath}: {e}")

        if not records:
            print(f"    ⚠️ Nenhum haplótipo encontrado para {gene} nos ficheiros processados.")
            continue

        # --- CORREÇÃO DO BUG 2 ---
        # O nome do ficheiro de saída deve ser "GENE.fasta", não "GENE_all_isolates.fasta"
        # para ser compatível com o próximo script (haplotype_analysis.py)
        output_path = os.path.join(results_dir, f"{gene}.fasta")
        SeqIO.write(records, output_path, "fasta")
        print(f"    ✅ {gene}: {len(records)} haplótipos gravados -> {output_path}")

# -------------------------------------------------------------
# Entry point
# -------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Uso: python haplotype_reconstruction.py [--reconstruct_single|--consolidate_all] <args> [--genes gene1,gene2,...]")

    # Parse optional --genes argument
    genes_arg = None
    if "--genes" in sys.argv:
        try:
            idx = sys.argv.index("--genes")
            genes_arg = sys.argv[idx + 1]
            sys.argv = sys.argv[:idx] + sys.argv[idx+2:] # remove --genes <arg>
        except IndexError:
            sys.exit("Erro: --genes requer um argumento (ex: ERG11,GSC1)")
        except ValueError:
            pass # --genes not found

    if sys.argv[1] == "--reconstruct_single" and len(sys.argv) >= 5:
        ref, resdir, isolate = sys.argv[2:5]
        genes = parse_gene_list(ref, genes_arg)
        reconstruct_single(ref, resdir, isolate, genes)

    elif sys.argv[1] == "--consolidate_all" and len(sys.argv) >= 4:
        # --- CORREÇÃO DO BUG 1 ---
        # Lemos o sys.argv[2] (results_dir) e sys.argv[3] (reference_fasta)
        resdir = sys.argv[2]
        fasta_path = sys.argv[3] # O caminho para o FASTA de referência

        if not os.path.exists(fasta_path):
            sys.exit(f"Erro: Ficheiro de referência FASTA não encontrado em {fasta_path}")

        # O parse_gene_list agora recebe o caminho correto
        genes = parse_gene_list(fasta_path, genes_arg)
        consolidate_all(resdir, genes)

    else:
        print("Erro: argumentos inválidos para o modo selecionado.")
        print("Modo --reconstruct_single requer: <reference_fasta> <results_dir> <isolate_name>")
        print("Modo --consolidate_all requer: <results_dir> <reference_fasta>")
        sys.exit(1)
