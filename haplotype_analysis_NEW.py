import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from tqdm import tqdm
from Bio.Align import PairwiseAligner, substitution_matrices
import argparse

# === CONFIGURATION AND SETUP ===

GENETIC_CODES = {
    1: "The Standard Code",
    2: "The Vertebrate Mitochondrial Code",
    3: "The Yeast Mitochondrial Code",
    4: "The Mold, Protozoan, and Coelenterate Mitochondrial Code and Mycoplasma/Spiroplasma Code",
    5: "The Invertebrate Mitochondrial Code",
    6: "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
    9: "The Echinoderm and Flatworm Mitochondrial Code",
    10: "The Euplotid Nuclear Code",
    11: "The Bacterial, Archaeal and Plant Plastid Code",
    12: "The Alternative Yeast Nuclear Code",
    13: "The Ascidian Mitochondrial Code",
    14: "The Alternative Flatworm Mitochondrial Code",
    15: "Blepharisma Nuclear Code",
    16: "Chlorophycean Mitochondrial Code",
    21: "Trematode Mitochondrial Code",
    22: "Scenedesmus obliquus Mitochondrial Code",
    23: "Thraustochytrium Mitochondrial Code",
    24: "Rhabdopleuridae Mitochondrial Code",
    25: "Candidate Division SR1 and Gracilibacteria Code",
    26: "Pachysolen tannophilus Nuclear Code",
    27: "Karyorelict Nuclear Code",
    28: "Condylostoma Nuclear Code",
    29: "Mesodinium Nuclear Code",
    30: "Peritrich Nuclear Code",
    31: "Blastocrithidia Nuclear Code",
    33: "Cephalodiscidae Mitochondrial UAA-Tyr Code",
}


def translate_with_trim(seq, table_id):
    """
    Translates a nucleotide sequence using Biopython's native fast translation.
    Truncates the protein sequence right after the first stop codon (*).
    """
    seq_str = str(seq).upper().replace("\n", "").replace(" ", "")
    # Ajusta para múltiplos de 3 para evitar avisos do Biopython
    remainder = len(seq_str) % 3
    if remainder != 0:
        seq_str = seq_str[:-remainder]

    if not seq_str:
        return ""

    bio_seq = Seq(seq_str)
    # Tradução nativa ultra rápida
    full_protein = str(bio_seq.translate(table=table_id))

    if "*" in full_protein:
        stop_index = full_protein.find("*")
        return full_protein[: stop_index + 1]  # Inclui o '*'

    return full_protein


def get_sort_key(variant_string):
    """
    Extracts the first integer position from a variant string for sorting.
    """
    match = re.match(r"^[A-Z\*](\d+)", variant_string, re.IGNORECASE)
    if match:
        return int(match.group(1))

    match = re.match(r"^[A-Z\*](\d+)_", variant_string, re.IGNORECASE)
    if match:
        return int(match.group(1))

    if variant_string.startswith("Frameshift"):
        match_fs = re.search(r"at position (\d+)", variant_string)
        if match_fs:
            return int(match_fs.group(1))

    return 99999


def get_variants_from_alignment(
    ref_nt, hap_nt, hap_id, nt_aligner, prot_aligner, table_id
):
    """
    Performs a robust two-step alignment using PairwiseAligner.
    """
    variants = set()

    # 1. Nucleotide Alignment (Global)
    nt_alignments = nt_aligner.align(ref_nt, hap_nt)
    if not nt_alignments:
        return variants, "", ""

    # Extrair o melhor alinhamento mapeado em formato de string
    alignment_info = nt_alignments[0]
    aln_ref_nt, aln_hap_nt = alignment_info[0], alignment_info[1]

    ref_nt_pos = 0
    in_indel = False
    frame_shift_accumulator = 0

    for r_char, h_char in zip(aln_ref_nt, aln_hap_nt):
        is_ref_gap = r_char == "-"
        is_hap_gap = h_char == "-"

        if is_ref_gap or is_hap_gap:
            if not in_indel:
                in_indel = True
                indel_start_pos = ref_nt_pos + 1
                indel_len = 0
                indel_type = "INS" if is_ref_gap else "DEL"

            if is_ref_gap:
                indel_len += 1
                frame_shift_accumulator = (frame_shift_accumulator + 1) % 3
            else:
                indel_len += 1
                frame_shift_accumulator = (frame_shift_accumulator - 1) % 3

        elif in_indel:
            in_indel = False
            if indel_len % 3 != 0:
                variants.add(
                    f"Frameshift {indel_type} of {indel_len} nt at position {indel_start_pos}"
                )

        if not is_ref_gap:
            ref_nt_pos += 1

    if frame_shift_accumulator != 0:
        variants.add(
            f"Catastrophic Frameshift: Frame not recovered (Net Shift: {frame_shift_accumulator} nt) in {hap_id}"
        )
        return (
            variants,
            translate_with_trim(ref_nt, table_id),
            translate_with_trim(hap_nt, table_id),
        )

    # 2. Protein Translation and Alignment
    ref_protein = translate_with_trim(ref_nt, table_id)
    hap_protein = translate_with_trim(hap_nt, table_id)

    if not ref_protein and not hap_protein:
        return variants, "", ""

    # Mudado para alinhamento global para evitar "comer" as pontas truncadas
    prot_alignments = prot_aligner.align(ref_protein, hap_protein)
    if not prot_alignments:
        return variants, ref_protein, hap_protein

    best_prot_aln = prot_alignments[0]
    aln_ref_prot, aln_hap_prot = best_prot_aln[0], best_prot_aln[1]

    aa_pos = 0
    i = 0
    aln_len = len(aln_ref_prot)

    while i < aln_len:
        r_char = aln_ref_prot[i]
        h_char = aln_hap_prot[i]

        if r_char != "-":
            aa_pos += 1

        if r_char != h_char:
            if r_char != "-" and h_char != "-":
                variants.add(f"{r_char}{aa_pos}{h_char}")
                i += 1

            elif r_char == "-":
                ins_seq = ""
                while i < aln_len and aln_ref_prot[i] == "-":
                    ins_seq += aln_hap_prot[i]
                    i += 1
                preceding_aa = ref_protein[aa_pos - 1] if aa_pos > 0 else "M"
                variants.add(f"{preceding_aa}{aa_pos}ins{ins_seq}")

            elif h_char == "-":
                del_seq = ""
                del_start_pos = aa_pos
                is_first = True
                while i < aln_len and aln_hap_prot[i] == "-":
                    del_seq += aln_ref_prot[i]
                    if not is_first:
                        aa_pos += 1
                    is_first = False
                    i += 1

                del_start_aa = del_seq[0] if del_seq else "?"

                if len(del_seq) == 1:
                    variants.add(f"{del_start_aa}{del_start_pos}del")
                else:
                    del_end_aa = del_seq[-1]
                    del_end_pos = del_start_pos + len(del_seq) - 1
                    variants.add(
                        f"{del_start_aa}{del_start_pos}_{del_end_aa}{del_end_pos}del"
                    )
        else:
            i += 1

    # Filtrar mutações após o primeiro stop codon do haplotipo
    first_stop_pos = -1
    if "*" in hap_protein:
        first_stop_pos = hap_protein.find("*") + 1

    if first_stop_pos != -1:
        variants_to_keep = set()
        for var in variants:
            if var.startswith("Frameshift") or var.startswith("Catastrophic"):
                variants_to_keep.add(var)
                continue

            var_pos = get_sort_key(var)
            if var_pos != 99999 and var_pos <= first_stop_pos:
                variants_to_keep.add(var)
        variants = variants_to_keep

    return variants, aln_ref_prot, aln_hap_prot


def extract_isolate_id(header):
    """
    Extracts isolate ID. Example: gene_Haplotype1_IsolateName -> IsolateName
    """
    first_word = header.split()[0]
    m = re.search(r"_(?:Haplotype[12])_(.+)$", first_word)
    if m:
        return m.group(1)
    return first_word


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze phased haplotype FASTA files against a reference sequence to identify protein-level variants.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "input_dir", type=str, help="Directory containing haplotype FASTA files"
    )
    parser.add_argument(
        "reference_file", type=str, help="Path to reference CDS FASTA file."
    )
    parser.add_argument(
        "output_dir", type=str, help="Directory where outputs will be saved."
    )
    parser.add_argument(
        "--transl_table",
        type=int,
        default=1,
        choices=GENETIC_CODES.keys(),
        help="Genetic code table number (default: 1).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    INPUT_DIR = args.input_dir
    REFERENCE_FILE = args.reference_file
    OUTPUT_DIR = args.output_dir
    transl_table_id = args.transl_table

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\n### Haplotype Mutation Analysis (PairwiseAligner Edition) ###")
    print(f" Input Directory: {INPUT_DIR}")
    print(f" Reference File: {REFERENCE_FILE}")
    print(f" Output Directory: {OUTPUT_DIR}\n")
    print(
        f" Translation Table: {transl_table_id} ({GENETIC_CODES[transl_table_id]})\n"
    )

    # Inicializar os novos Alinhadores modernos do Biopython 1.81+
    nt_aligner = PairwiseAligner()
    nt_aligner.mode = "global"
    nt_aligner.match_score = 2
    nt_aligner.mismatch_score = -1
    nt_aligner.open_gap_score = -10
    nt_aligner.extend_gap_score = -1

    prot_aligner = PairwiseAligner()
    prot_aligner.mode = "global"  # Mudado para global para apanhar pontas truncadas
    prot_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    prot_aligner.open_gap_score = -10
    prot_aligner.extend_gap_score = -1

    try:
        ref_records = SeqIO.to_dict(SeqIO.parse(REFERENCE_FILE, "fasta"))
    except FileNotFoundError:
        print(f"ERROR: Reference file not found at {REFERENCE_FILE}")
        return

    for fasta_file in tqdm(sorted(os.listdir(INPUT_DIR))):
        full_path = os.path.join(INPUT_DIR, fasta_file)
        if not os.path.isfile(full_path):
            continue

        if not fasta_file.endswith(".fasta") or fasta_file == os.path.basename(
            REFERENCE_FILE
        ):
            continue

        gene_name = fasta_file.replace(".fasta", "")

        ref_rec = ref_records.get(gene_name)
        if ref_rec is None:
            print(
                f" Warning: Reference not found for gene {gene_name} in {REFERENCE_FILE}, skipping."
            )
            continue

        ref_seq = str(ref_rec.seq).upper()
        print(f"\nProcessing {gene_name}...")

        try:
            records = list(SeqIO.parse(full_path, "fasta"))
        except Exception as e:
            print(
                f" Could not read FASTA file: {fasta_file}. Skipping. Error: {e}"
            )
            continue

        isolates = {}
        anomalies = []

        for rec in records:
            isolate = extract_isolate_id(rec.id)
            isolates.setdefault(isolate, []).append(str(rec.seq).upper())

        all_variants_in_gene = set()
        rows = []
        protein_alignments_log = []

        for isolate, seqs in isolates.items():
            if len(seqs) != 2:
                anomalies.append(
                    f"{isolate}: Found {len(seqs)} haplotypes (expected 2)"
                )
                continue

            vars_h1, aln_ref_h1, aln_hap_h1 = get_variants_from_alignment(
                ref_seq,
                seqs[0],
                "Haplotype1",
                nt_aligner,
                prot_aligner,
                transl_table_id,
            )
            vars_h2, aln_ref_h2, aln_hap_h2 = get_variants_from_alignment(
                ref_seq,
                seqs[1],
                "Haplotype2",
                nt_aligner,
                prot_aligner,
                transl_table_id,
            )

            all_vars_for_isolate = vars_h1.union(vars_h2)
            all_variants_in_gene.update(all_vars_for_isolate)

            isolate_data = {"Isolate": isolate}
            for var in all_vars_for_isolate:
                in_h1 = var in vars_h1
                in_h2 = var in vars_h2
                isolate_data[var] = "Hom" if in_h1 and in_h2 else "Het"

            rows.append(isolate_data)

            if all_vars_for_isolate:
                if vars_h1:
                    protein_alignments_log.append(
                        (
                            isolate,
                            "Haplotype1",
                            aln_ref_h1,
                            aln_hap_h1,
                            sorted(list(vars_h1), key=get_sort_key),
                        )
                    )
                if vars_h2:
                    protein_alignments_log.append(
                        (
                            isolate,
                            "Haplotype2",
                            aln_ref_h2,
                            aln_hap_h2,
                            sorted(list(vars_h2), key=get_sort_key),
                        )
                    )

        if not rows:
            df = pd.DataFrame(columns=["Isolate"])
        else:
            df = pd.DataFrame(rows).fillna("")

        all_variants_list = list(all_variants_in_gene)
        sorted_variants = sorted(all_variants_list, key=get_sort_key)
        cols = ["Isolate"] + sorted_variants
        df = df.reindex(columns=cols)

        out_csv = os.path.join(OUTPUT_DIR, f"{gene_name}_variants_table.csv")
        out_report = os.path.join(OUTPUT_DIR, f"{gene_name}_variant_summary.txt")
        out_align = os.path.join(
            OUTPUT_DIR, f"{gene_name}_protein_alignments.txt"
        )

        df.to_csv(out_csv, index=False, sep=";")

        with open(out_align, "w", encoding="utf-8") as f_align:
            f_align.write(f"### Detailed Protein Alignments for {gene_name}\n")
            f_align.write(
                "NOTE: Use this file to debug indel/substitution ambiguity.\n\n"
            )

            for (
                isolate,
                hap_name,
                aln_ref,
                aln_hap,
                variants_list,
            ) in protein_alignments_log:
                f_align.write("=" * 80 + "\n")
                f_align.write(f"Isolate: {isolate} | Haplotype: {hap_name}\n")
                f_align.write(f"Detected Variants: {', '.join(variants_list)}\n\n")

                block_size = 70
                for start in range(0, len(aln_ref), block_size):
                    end = start + block_size
                    f_align.write(f"Ref: {aln_ref[start:end]}\n")
                    f_align.write(f"Hap: {aln_hap[start:end]}\n")
                    f_align.write("\n")
                f_align.write("\n")

        with open(out_report, "w", encoding="utf-8") as f:
            f.write(f"### Variant Summary for {gene_name}\n\n")
            f.write(f"Total isolates processed: {len(isolates)}\n")
            f.write(
                f"Total unique variants detected: {len(all_variants_in_gene)}\n\n"
            )

            aa_variants = [
                v
                for v in sorted_variants
                if not v.startswith("Frameshift")
                and not v.startswith("Catastrophic")
            ]
            nt_events = [
                v
                for v in sorted_variants
                if v.startswith("Frameshift") or v.startswith("Catastrophic")
            ]

            f.write(
                "Detected AA variants (Substitutions, In-frame Indels, Stops):\n"
            )
            for v in aa_variants:
                f.write(f" - {v}\n")
            f.write("\nDetected NT (Frameshift/Catastrophic) variants:\n")
            for v in nt_events:
                f.write(f" - {v}\n")
            f.write("\nAnomalies and warnings:\n")
            if anomalies:
                for a in anomalies:
                    f.write(f" Alert: {a}\n")
            else:
                f.write(" None detected.\n")

        print(
            f"✅ {gene_name} done! ({len(isolates)} isolates, {len(all_variants_in_gene)} variants)"
        )

    print("\n Analysis complete! All results saved in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
import osimport reimport pandas as pdfrom Bio import SeqIOfrom Bio.Data import CodonTablefrom Bio.Seq import Seqfrom tqdm import tqdmfrom Bio.Align import PairwiseAligner, substitution_matricesimport argparse=== CONFIGURATION AND SETUP ===Genetic codes referenceGENETIC_CODES = {1: "The Standard Code",2: "The Vertebrate Mitochondrial Code",3: "The Yeast Mitochondrial Code",4: "The Mold, Protozoan, and Coelenterate Mitochondrial Code and Mycoplasma/Spiroplasma Code",5: "The Invertebrate Mitochondrial Code",6: "The Ciliate, Dasycladacean and Hexamita Nuclear Code",9: "The Echinoderm and Flatworm Mitochondrial Code",10: "The Euplotid Nuclear Code",11: "The Bacterial, Archaeal and Plant Plastid Code",12: "The Alternative Yeast Nuclear Code",13: "The Ascidian Mitochondrial Code",14: "The Alternative Flatworm Mitochondrial Code",15: "Blepharisma Nuclear Code",16: "Chlorophycean Mitochondrial Code",21: "Trematode Mitochondrial Code",22: "Scenedesmus obliquus Mitochondrial Code",23: "Thraustochytrium Mitochondrial Code",24: "Rhabdopleuridae Mitochondrial Code",25: "Candidate Division SR1 and Gracilibacteria Code",26: "Pachysolen tannophilus Nuclear Code",27: "Karyorelict Nuclear Code",28: "Condylostoma Nuclear Code",29: "Mesodinium Nuclear Code",30: "Peritrich Nuclear Code",31: "Blastocrithidia Nuclear Code",33: "Cephalodiscidae Mitochondrial UAA-Tyr Code",}def translate_with_trim(seq, transl_table_id):"""Translate a nucleotide sequence using the selected genetic code.Stops translation at the first internal stop codon, returning '*' andtruncating the resulting protein sequence at that point.Optimized to use Biopython's native fast C-implementation."""seq_str = str(seq).upper().replace("\n", "").replace(" ", "")# Use native Biopython translation (much faster and robust)
protein = str(Seq(seq_str).translate(table=transl_table_id))

if '*' in protein:
    # Truncate at first stop but keep the '*' as in the original logic
    protein = protein.split('*')[0] + '*'
    
return protein
def get_sort_key(variant_string):"""Helper function to extract the first integer (position) from a variant stringfor numerical sorting of output columns (e.g., G159V -> 159)."""match = re.match(r'^A-Z*', variant_string, re.IGNORECASE)if match:return int(match.group(1))match = re.match(r'^[A-Z\*](\d+)_', variant_string, re.IGNORECASE)
if match:
    return int(match.group(1))

if variant_string.startswith("Frameshift"):
    match_fs = re.search(r'at position (\d+)', variant_string)
    if match_fs:
        return int(match_fs.group(1))

return 99999
def get_variants_from_alignment(ref_nt, hap_nt, hap_id, transl_table_id, aligner_nt, aligner_prot):"""Performs a robust two-step alignment (NT then Protein) to identify variants.1. Nucleotide alignment checks for frameshifts.2. Protein alignment checks for substitutions and in-frame indels.3. Filters out all AA changes that occur after the first stop codon in the haplotype.Returns: (variants, aligned_ref_protein, aligned_hap_protein)"""variants = set()# Modern Biopython pairwise alignment using PairwiseAligner
nt_alignments = aligner_nt.align(ref_nt, hap_nt)
if not nt_alignments:
    return variants, "", ""

# Unpack the first alignment with gaps using standard Align slicing
aln_ref_nt = str(nt_alignments[0][0, :])
aln_hap_nt = str(nt_alignments[0][1, :])

ref_nt_pos = 0
in_indel = False
frame_shift_accumulator = 0

for r_char, h_char in zip(aln_ref_nt, aln_hap_nt):
    is_ref_gap = (r_char == '-')
    is_hap_gap = (h_char == '-')

    if is_ref_gap or is_hap_gap:
        if not in_indel:
            in_indel = True
            indel_start_pos = ref_nt_pos + 1
            indel_len = 0
            indel_type = "INS" if is_ref_gap else "DEL"

        if is_ref_gap:
            indel_len += 1
            frame_shift_accumulator = (frame_shift_accumulator + 1) % 3
        else:
            indel_len += 1
            frame_shift_accumulator = (frame_shift_accumulator - 1) % 3

    elif in_indel:
        in_indel = False
        if indel_len % 3 != 0:
            variants.add(
                f"Frameshift {indel_type} of {indel_len} nt at position {indel_start_pos}"
            )

    if not is_ref_gap:
        ref_nt_pos += 1

if frame_shift_accumulator != 0:
    variants.add(
        f"Catastrophic Frameshift: Frame not recovered (Net Shift: {frame_shift_accumulator} nt) in {hap_id}"
    )
    return variants, translate_with_trim(ref_nt, transl_table_id), translate_with_trim(hap_nt, transl_table_id)

ref_protein = translate_with_trim(ref_nt, transl_table_id)
hap_protein = translate_with_trim(hap_nt, transl_table_id)

if not ref_protein and not hap_protein:
    return variants, "", ""

# Align protein sequences using modern PairwiseAligner (Global mode ensures no edge mutations are clipped)
prot_alignments = aligner_prot.align(ref_protein, hap_protein)
if not prot_alignments:
    return variants, ref_protein, hap_protein

# Get gapped sequence strings
aln_ref_prot = str(prot_alignments[0][0, :])
aln_hap_prot = str(prot_alignments[0][1, :])

aa_pos = 0
i = 0
aln_len = len(aln_ref_prot)

while i < aln_len:
    r_char = aln_ref_prot[i]
    h_char = aln_hap_prot[i]

    if r_char != '-':
        aa_pos += 1

    if r_char != h_char:
        if r_char != '-' and h_char != '-':
            variants.add(f"{r_char}{aa_pos}{h_char}")
            i += 1

        elif r_char == '-':
            ins_seq = ""
            while i < aln_len and aln_ref_prot[i] == '-':
                ins_seq += aln_hap_prot[i]
                i += 1
            preceding_aa = ref_protein[aa_pos - 1] if aa_pos > 0 else 'M'
            variants.add(f"{preceding_aa}{aa_pos}ins{ins_seq}")

        elif h_char == '-':
            del_seq = ""
            del_start_pos = aa_pos
            is_first = True
            while i < aln_len and aln_hap_prot[i] == '-':
                del_seq += aln_ref_prot[i]
                if not is_first:
                     aa_pos += 1
                is_first = False
                i += 1

            del_start_aa = del_seq[0] if del_seq else '?'

            if len(del_seq) == 1:
                variants.add(f"{del_start_aa}{del_start_pos}del")
            else:
                del_end_aa = del_seq[-1]
                del_end_pos = del_start_pos + len(del_seq) - 1
                variants.add(f"{del_start_aa}{del_start_pos}_{del_end_aa}{del_end_pos}del")
    else:
        i += 1

first_stop_pos = -1
if '*' in hap_protein:
    first_stop_pos = hap_protein.find('*') + 1

if first_stop_pos != -1:
    variants_to_keep = set()
    for var in variants:
        if var.startswith("Frameshift") or var.startswith("Catastrophic"):
             variants_to_keep.add(var)
             continue

        var_pos = get_sort_key(var)

        if var_pos != 99999 and var_pos <= first_stop_pos:
            variants_to_keep.add(var)
    variants = variants_to_keep

return variants, aln_ref_prot, aln_hap_prot
def extract_isolate_id(header):"""Extracts the isolate ID from a header of type:Haplotype1<isolate_id>or Haplotype2<isolate_id>Keeps underscores present in the isolate_id (e.g., A_19_3)."""first_word = header.split()[0]m = re.search(r'(?:Haplotype[12])(.+)$', first_word)if m:return m.group(1)return first_worddef parse_args():parser = argparse.ArgumentParser(description="Analyze phased haplotype FASTA files against a reference sequence to identify protein-level variants (substitutions, indels, and frameshifts).",formatter_class=argparse.RawTextHelpFormatter)parser.add_argument('input_dir', type=str, help="Directory containing the haplotype FASTA files")parser.add_argument('reference_file', type=str, help="Path to the single FASTA file containing all reference CDS sequences.")parser.add_argument('output_dir', type=str, help="Directory where the output CSV tables and text reports will be saved.")parser.add_argument('--transl_table', type=int, default=1, choices=GENETIC_CODES.keys(), help="Genetic code table number (default: 1).")return parser.parse_args()=== MAIN PIPELINE ===def main():args = parse_args()INPUT_DIR = args.input_dirREFERENCE_FILE = args.reference_fileOUTPUT_DIR = args.output_dirtransl_table_id = args.transl_tableos.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n### Haplotype Mutation Analysis (Variant Caller) ###")
print(f" Input Directory: {INPUT_DIR}")
print(f" Reference File: {REFERENCE_FILE}")
print(f" Output Directory: {OUTPUT_DIR}\n")

print(f" Translation Table: {transl_table_id} ({GENETIC_CODES[transl_table_id]})\n")

# Load BLOSUM62 once
try:
    blosum62 = substitution_matrices.load("BLOSUM62")
except Exception as e:
    print(f"ERROR: Could not load BLOSUM62 matrix. {e}")
    return

# Define Aligners once to maximize speed and efficiency
aligner_nt = PairwiseAligner()
aligner_nt.mode = "global"
aligner_nt.match_score = 2
aligner_nt.mismatch_score = -1
aligner_nt.open_gap_score = -10
aligner_nt.extend_gap_score = -1

aligner_prot = PairwiseAligner()
aligner_prot.mode = "global"  # Corrected to global to prevent clipping of variants at termini
aligner_prot.substitution_matrix = blosum62
aligner_prot.open_gap_score = -10
aligner_prot.extend_gap_score = -1

try:
    ref_records = SeqIO.to_dict(SeqIO.parse(REFERENCE_FILE, "fasta"))
except FileNotFoundError:
    print(f"ERROR: Reference file not found at {REFERENCE_FILE}")
    return

for fasta_file in tqdm(sorted(os.listdir(INPUT_DIR))):
    full_path = os.path.join(INPUT_DIR, fasta_file)
    if not os.path.isfile(full_path):
        continue

    if not fasta_file.endswith(".fasta") or fasta_file == os.path.basename(REFERENCE_FILE):
        continue

    gene_name = fasta_file.replace(".fasta", "")

    ref_rec = ref_records.get(gene_name)
    if ref_rec is None:
        ref_seq = ""
    else:
        ref_seq = str(ref_rec.seq).upper()

    if not ref_seq:
        print(f" Warning: Reference not found for gene **{gene_name}** in {REFERENCE_FILE}, skipping.")
        continue

    print(f"\nProcessing {gene_name}...")

    try:
        records = list(SeqIO.parse(full_path, "fasta"))
    except Exception as e:
        print(f" Could not read FASTA file: {fasta_file}. Skipping. Error: {e}")
        continue

    isolates = {}
    anomalies = []

    for rec in records:
        isolate = extract_isolate_id(rec.id)
        isolates.setdefault(isolate, []).append(str(rec.seq).upper())

    all_variants_in_gene = set()
    rows = []
    protein_alignments_log = []

    for isolate, seqs in isolates.items():
        if len(seqs) != 2:
            anomalies.append(f"{isolate}: Found {len(seqs)} haplotypes (expected 2)")
            continue

        vars_h1, aln_ref_h1, aln_hap_h1 = get_variants_from_alignment(
            ref_seq, seqs[0], "Haplotype1", transl_table_id, aligner_nt, aligner_prot
        )
        vars_h2, aln_ref_h2, aln_hap_h2 = get_variants_from_alignment(
            ref_seq, seqs[1], "Haplotype2", transl_table_id, aligner_nt, aligner_prot
        )

        all_vars_for_isolate = vars_h1.union(vars_h2)
        all_variants_in_gene.update(all_vars_for_isolate)

        isolate_data = {"Isolate": isolate}

        for var in all_vars_for_isolate:
            in_h1 = var in vars_h1
            in_h2 = var in vars_h2
            isolate_data[var] = "Hom" if in_h1 and in_h2 else "Het"

        rows.append(isolate_data)

        if all_vars_for_isolate:
            if vars_h1:
                protein_alignments_log.append(
                    (isolate, "Haplotype1", aln_ref_h1, aln_hap_h1, sorted(list(vars_h1), key=get_sort_key))
                )
            if vars_h2:
                protein_alignments_log.append(
                    (isolate, "Haplotype2", aln_ref_h2, aln_hap_h2, sorted(list(vars_h2), key=get_sort_key))
                )

    if rows:
        df = pd.DataFrame(rows).fillna("")
        all_variants_list = list(all_variants_in_gene)
        sorted_variants = sorted(all_variants_list, key=get_sort_key)
        cols = ["Isolate"] + sorted_variants
        df = df.reindex(columns=cols, fill_value="")
    else:
        df = pd.DataFrame(columns=["Isolate"])

    out_csv = os.path.join(OUTPUT_DIR, f"{gene_name}_variants_table.csv")
    out_report = os.path.join(OUTPUT_DIR, f"{gene_name}_variant_summary.txt")
    out_align = os.path.join(OUTPUT_DIR, f"{gene_name}_protein_alignments.txt")

    df.to_csv(out_csv, index=False, sep=';')

    with open(out_align, "w", encoding="utf-8") as f_align:
        f_align.write(f"### Detailed Protein Alignments for {gene_name}\n")
        f_align.write("NOTE: Use this file to debug indel/substitution ambiguity, such as positions after a deletion.\n\n")

        for isolate, hap_name, aln_ref, aln_hap, variants_list in protein_alignments_log:
            f_align.write("=" * 80 + "\n")
            f_align.write(f"Isolate: {isolate} | Haplotype: {hap_name}\n")
            f_align.write(f"Detected Variants: {', '.join(variants_list)}\n\n")

            block_size = 70
            for start in range(0, len(aln_ref), block_size):
                end = start + block_size
                f_align.write(f"Ref: {aln_ref[start:end]}\n")
                f_align.write(f"Hap: {aln_hap[start:end]}\n")
                f_align.write("\n")
            f_align.write("\n")

    with open(out_report, "w", encoding="utf-8") as f:
        f.write(f"### Variant Summary for {gene_name}\n\n")
        f.write(f"Total isolates processed: {len(isolates)}\n")
        f.write(f"Total unique variants detected: {len(all_variants_in_gene)}\n\n")

        sorted_variants = sorted(list(all_variants_in_gene), key=get_sort_key)
        aa_variants = [v for v in sorted_variants if not v.startswith("Frameshift") and not v.startswith("Catastrophic")]
        nt_events = [v for v in sorted_variants if v.startswith("Frameshift") or v.startswith("Catastrophic")]

        f.write("Detected AA variants (Substitutions, In-frame Indels, Stops):\n")
        for v in aa_variants:
            f.write(f" - {v}\n")
        f.write("\nDetected NT (Frameshift/Catastrophic) variants:\n")
        for v in nt_events:
            f.write(f" - {v}\n")
        f.write("\nAnomalies and warnings:\n")
        if anomalies:
            for a in anomalies:
                f.write(f" Alert: {a}\n")
        else:
            f.write(" None detected.\n")

    print(f"✅ {gene_name} done! ({len(isolates)} isolates, {len(all_variants_in_gene)} variants)")
    print(f"   Saved table to {os.path.basename(out_csv)} and alignment log to {os.path.basename(out_align)}")

print("\n Analysis complete! All results saved in:", OUTPUT_DIR)
if name == "main":main()
