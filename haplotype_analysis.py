import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from tqdm import tqdm
from Bio import pairwise2
from Bio.Align import substitution_matrices
import argparse

# === CONFIGURATION AND SETUP ===

# Load BLOSUM62 matrix using the modern Biopython module
BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Translation table 1 — Standard Nuclear Code
translation_table = CodonTable.unambiguous_dna_by_id[1]

# === HELPER FUNCTIONS ===
def translate_with_trim(seq):
    """
    Translate a nucleotide sequence using yeast code (table 12).
    Stops translation at the first internal stop codon, returning '*' and
    truncating the resulting protein sequence at that point.
    """
    protein = []
    # Ensure seq is a string for .upper()
    seq_str = str(seq).upper().replace("\n", "").replace(" ", "")

    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i + 3]
        if len(codon) < 3:
            break

        aa = translation_table.forward_table.get(codon, 'X')

        if codon in translation_table.stop_codons:
            protein.append('*')
            break  # Truncate at stop

        protein.append(aa)
    return "".join(protein)


def get_sort_key(variant_string):
    """
    Helper function to extract the first integer (position) from a variant string
    for numerical sorting of output columns (e.g., G159V -> 159).
    """
    # Standard AA change: G159V or G159*
    match = re.match(r'^[A-Z\*](\d+)', variant_string, re.IGNORECASE)
    if match:
        return int(match.group(1))

    # In-frame indel: G85_A86insL or G85_K87del
    match = re.match(r'^[A-Z\*](\d+)_', variant_string, re.IGNORECASE)
    if match:
        return int(match.group(1))

    # Frameshift: Frameshift DEL of 14 nt at position 420...
    if variant_string.startswith("Frameshift"):
        match_fs = re.search(r'at position (\d+)', variant_string)
        if match_fs:
            return int(match_fs.group(1))

    # Return a high number for un-numbered variants so they go last
    return 99999


def get_variants_from_alignment(ref_nt, hap_nt, hap_id):
    """
    Performs a robust two-step alignment (NT then Protein) to identify variants.

    1. Nucleotide alignment checks for frameshifts.
    2. Protein alignment checks for substitutions and in-frame indels.
    3. Filters out all AA changes that occur after the first stop codon in the haplotype.

    Returns: (variants, aligned_ref_protein, aligned_hap_protein)
    """
    variants = set()

    # --- Step 1: Nucleotide Alignment (Frameshift Check with Recovery Logic) ---
    nt_alignments = pairwise2.align.globalms(ref_nt, hap_nt, 2, -1, -10, -1)
    if not nt_alignments:
        return variants, "", ""

    aln_ref_nt, aln_hap_nt = nt_alignments[0].seqA, nt_alignments[0].seqB

    ref_nt_pos = 0
    in_indel = False
    # Tracks the cumulative frame shift (0, 1, or 2). Must be 0 at the end.
    frame_shift_accumulator = 0

    # Scan the NT alignment for any indel not divisible by 3
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
            # Check if this single block caused a non-in-frame change
            if indel_len % 3 != 0:
                variants.add(
                    f"Frameshift {indel_type} of {indel_len} nt at position {indel_start_pos}"
                )

        if not is_ref_gap:
            ref_nt_pos += 1

    # Final check for catastrophic frameshift after all indel blocks
    if frame_shift_accumulator != 0:
        variants.add(
            f"Catastrophic Frameshift: Frame not recovered (Net Shift: {frame_shift_accumulator} nt) in {hap_id}"
        )
        # Catastrophic event, stop here.
        return variants, translate_with_trim(ref_nt), translate_with_trim(hap_nt)

    # --- Step 2: Protein Alignment (Substitutions & In-Frame Indels) ---

    ref_protein = translate_with_trim(ref_nt)
    hap_protein = translate_with_trim(hap_nt)

    if not ref_protein and not hap_protein:
        return variants, "", ""

    # Use local alignment with the BLOSUM62 substitution matrix.
    prot_alignments = pairwise2.align.localds(ref_protein, hap_protein, BLOSUM62, -10, -1)
    if not prot_alignments:
        return variants, ref_protein, hap_protein

    aln_ref_prot, aln_hap_prot = prot_alignments[0].seqA, prot_alignments[0].seqB

    aa_pos = 0 # Tracks position in the REFERENCE sequence
    i = 0      # Tracks index in the ALIGNED sequence
    aln_len = len(aln_ref_prot)

    while i < aln_len:
        r_char = aln_ref_prot[i]
        h_char = aln_hap_prot[i]

        if r_char != '-':
            aa_pos += 1 # Advance reference position if current aligned position is not a gap in Ref

        if r_char != h_char:

            if r_char != '-' and h_char != '-':
                # Case 1: Substitution (Mismatched but not a gap)
                variants.add(f"{r_char}{aa_pos}{h_char}")
                i += 1

            elif r_char == '-':
                # Case 2: Insertion (Gap in reference)
                ins_seq = ""
                while i < aln_len and aln_ref_prot[i] == '-':
                    ins_seq += aln_hap_prot[i]
                    i += 1

                # Insertion reporting is relative to the *preceding* amino acid (aa_pos)
                preceding_aa = ref_protein[aa_pos - 1] if aa_pos > 0 else 'M'
                variants.add(f"{preceding_aa}{aa_pos}ins{ins_seq}")

            elif h_char == '-':
                # Case 3: Deletion (Gap in haplotype)
                del_seq = ""
                del_start_pos = aa_pos

                is_first = True
                while i < aln_len and aln_hap_prot[i] == '-':
                    del_seq += aln_ref_prot[i]
                    if not is_first:
                         aa_pos += 1
                    is_first = False
                    i += 1

                del_start_aa = del_seq[0]

                if len(del_seq) == 1:
                    # Single AA deletion (e.g., L144del)
                    variants.add(f"{del_start_aa}{del_start_pos}del")
                else:
                    # Multi AA deletion (e.g., L144_T145del)
                    del_end_aa = del_seq[-1]
                    del_end_pos = del_start_pos + len(del_seq) - 1
                    variants.add(f"{del_start_aa}{del_start_pos}_{del_end_aa}{del_end_pos}del")
        else:
            # Case 4: Match
            i += 1

    # --- Logic: Filter variants after a stop codon in the HAPLOTYPE ---
    first_stop_pos = -1
    if '*' in hap_protein:
        first_stop_pos = hap_protein.find('*') + 1

    if first_stop_pos != -1:
        variants_to_keep = set()
        for var in variants:
            # Keep frameshift variants as they describe the cause of the truncation
            if var.startswith("Frameshift") or var.startswith("Catastrophic"):
                 variants_to_keep.add(var)
                 continue

            var_pos = get_sort_key(var)

            # Keep the variant if it's before or exactly at the stop position
            if var_pos != 99999 and var_pos <= first_stop_pos:
                variants_to_keep.add(var)
        variants = variants_to_keep

    return variants, aln_ref_prot, aln_hap_prot


def extract_isolate_id(header):
    """
    Extract isolate ID from FASTA header.
    (e.g., GSC1_Haplotype1_A_19_13 Alelo 1... -> A_19_13)
    """
    first_word = header.split()[0]
    parts = re.split(r'_Haplotype[12]_', first_word, maxsplit=1)
    if len(parts) > 1:
        # Returns the part after the haplotype indicator
        return parts[1].split("_")[0] # Assumes the isolate ID is the first part after the split
    return "Unknown"


def parse_args():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze phased haplotype FASTA files against a reference sequence to identify protein-level variants (substitutions, indels, and frameshifts).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'input_dir',
        type=str,
        help="Directory containing the haplotype FASTA files (e.g., 'ERG11.fasta', 'GSC1.fasta')."
    )
    parser.add_argument(
        'reference_file',
        type=str,
        help="Path to the single FASTA file containing all reference CDS sequences."
    )
    parser.add_argument(
        'output_dir',
        type=str,
        help="Directory where the output CSV tables and text reports will be saved."
    )
    return parser.parse_args()


# === MAIN PIPELINE ===
def main():
    args = parse_args()
    INPUT_DIR = args.input_dir
    REFERENCE_FILE = args.reference_file
    OUTPUT_DIR = args.output_dir

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\n### Haplotype Mutation Analysis (Variant Caller) ###")
    print(f" Input Directory: {INPUT_DIR}")
    print(f" Reference File: {REFERENCE_FILE}")
    print(f" Output Directory: {OUTPUT_DIR}\n")

    # Load reference sequences
    try:
        ref_records = SeqIO.to_dict(SeqIO.parse(REFERENCE_FILE, "fasta"))
    except FileNotFoundError:
        print(f"ERROR: Reference file not found at {REFERENCE_FILE}")
        return

    # Iterate over each gene FASTA in the input directory
    for fasta_file in tqdm(sorted(os.listdir(INPUT_DIR))):
        full_path = os.path.join(INPUT_DIR, fasta_file)

        # IMPROVEMENT: Skip if it is a directory or if it is not a file (symlink, etc.)
        if not os.path.isfile(full_path):
            continue

        # Check 1: Must be a .fasta file AND not the reference file
        if not fasta_file.endswith(".fasta") or fasta_file == os.path.basename(REFERENCE_FILE):
            continue

        gene_name = fasta_file.replace(".fasta", "")
        ref_seq = str(ref_records.get(gene_name, Seq("")).seq).upper()

        if not ref_seq:
            print(f" Warning: Reference not found for gene **{gene_name}** in {REFERENCE_FILE}, skipping.")
            continue

        print(f"\nProcessing {gene_name}...")

        try:
            records = list(SeqIO.parse(full_path, "fasta"))
        except:
            print(f" Could not read FASTA file: {fasta_file}. Skipping.")
            continue

        isolates = {}
        anomalies = []

        # Group haplotypes per isolate
        for rec in records:
            isolate = extract_isolate_id(rec.id)
            isolates.setdefault(isolate, []).append(str(rec.seq).upper())

        # Data storage
        all_variants_in_gene = set()
        rows = []
        protein_alignments_log = []

        for isolate, seqs in isolates.items():
            if len(seqs) != 2:
                anomalies.append(f"{isolate}: Found {len(seqs)} haplotypes (expected 2)")
                continue

            # --- Compare each haplotype to ref individually ---
            vars_h1, aln_ref_h1, aln_hap_h1 = get_variants_from_alignment(ref_seq, seqs[0], "Haplotype1")
            vars_h2, aln_ref_h2, aln_hap_h2 = get_variants_from_alignment(ref_seq, seqs[1], "Haplotype2")

            all_vars_for_isolate = vars_h1.union(vars_h2)
            all_variants_in_gene.update(all_vars_for_isolate)

            isolate_data = {"Isolate": isolate}

            # Determine zygosity (Het/Hom)
            for var in all_vars_for_isolate:
                in_h1 = var in vars_h1
                in_h2 = var in vars_h2

                # Het: In one haplotype. Hom: In both haplotypes.
                isolate_data[var] = "Hom" if in_h1 and in_h2 else "Het"

            rows.append(isolate_data)

            # --- Logging alignments if variants are found ---
            if all_vars_for_isolate:
                if vars_h1:
                    protein_alignments_log.append(
                        (isolate, "Haplotype1", aln_ref_h1, aln_hap_h1, sorted(list(vars_h1), key=get_sort_key))
                    )
                if vars_h2:
                    protein_alignments_log.append(
                        (isolate, "Haplotype2", aln_ref_h2, aln_hap_h2, sorted(list(vars_h2), key=get_sort_key))
                    )


        # Build DataFrame
        df = pd.DataFrame(rows).fillna("")

        # Sort columns: Isolate first, then all variants sorted by position
        all_variants_list = list(all_variants_in_gene)
        sorted_variants = sorted(all_variants_list, key=get_sort_key)

        cols = ["Isolate"] + sorted_variants
        df = df.reindex(columns=cols)

        # Save outputs
        out_csv = os.path.join(OUTPUT_DIR, f"{gene_name}_variants_table.csv")
        out_report = os.path.join(OUTPUT_DIR, f"{gene_name}_variant_summary.txt")
        out_align = os.path.join(OUTPUT_DIR, f"{gene_name}_protein_alignments.txt")

        df.to_csv(out_csv, index=False, sep=';') # Use semicolon for Excel/spreadsheet compatibility

        # --- Alignment Log Output ---
        with open(out_align, "w", encoding="utf-8") as f_align:
            f_align.write(f"### Detailed Protein Alignments for {gene_name}\n")
            f_align.write("NOTE: Use this file to debug indel/substitution ambiguity, such as positions after a deletion.\n\n")

            for isolate, hap_name, aln_ref, aln_hap, variants_list in protein_alignments_log:
                f_align.write("=" * 80 + "\n")
                f_align.write(f"Isolate: {isolate} | Haplotype: {hap_name}\n")
                f_align.write(f"Detected Variants: {', '.join(variants_list)}\n\n")

                # Print alignment in blocks for readability
                block_size = 70
                for start in range(0, len(aln_ref), block_size):
                    end = start + block_size
                    f_align.write(f"Ref: {aln_ref[start:end]}\n")
                    f_align.write(f"Hap: {aln_hap[start:end]}\n")
                    f_align.write("\n")

                f_align.write("\n")

        # --- Summary Report Output ---
        with open(out_report, "w", encoding="utf-8") as f:
            f.write(f"### Variant Summary for {gene_name}\n\n")
            f.write(f"Total isolates processed: {len(isolates)}\n")
            f.write(f"Total unique variants detected: {len(all_variants_in_gene)}\n\n")

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


if __name__ == "__main__":
    main()
