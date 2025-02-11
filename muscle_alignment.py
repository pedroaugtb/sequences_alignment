#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import zipfile

# Configuration
DATA_DIR = Path("ncbi_dataset/ncbi_dataset/data")
RESULTS_DIR = Path("results_muscle")
COLLECTED_FASTA = RESULTS_DIR / "all_sequences.faa"
ALIGNED_FASTA = RESULTS_DIR / "aligned.faa"
CONSERVE_ARCHIVE = RESULTS_DIR / "conserved_sequences.zip"
NON_CONSERVE_ARCHIVE = RESULTS_DIR / "non_conserved_sequences.zip"
AMALGAM_FASTA = RESULTS_DIR / "amalgam_sequence.faa"
AMALGAM_ARCHIVE = RESULTS_DIR / "amalgam_sequence.zip"
CONSERVATION_THRESHOLD = 0.11  # Increased threshold for stricter conservation

def collect_sequences(data_dir: Path, output_fasta: Path):
    """
    Collects all protein sequences from the data directory and writes them to a single FASTA file.
    Accumulates warnings about missing files and prints a summary at the end.
    """
    sequences = []
    missing_files = []  # List to store missing files

    for genome_dir in data_dir.iterdir():
        if genome_dir.is_dir():
            protein_faa = genome_dir / "protein.faa"
            if protein_faa.exists():
                for record in SeqIO.parse(protein_faa, "fasta"):
                    sequences.append(record)
            else:
                missing_files.append(protein_faa)

    if not sequences:
        print("Error: No sequences found.", file=sys.stderr)
        sys.exit(1)

    # Ensure the output directory exists
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"Collected {len(sequences)} sequences in {output_fasta}")

    # Summary of missing files
    if missing_files:
        print(f"Warning: {len(missing_files)} 'protein.faa' files were not found.", file=sys.stderr)
    else:
        print("All 'protein.faa' files were found.")

def run_muscle(input_fasta: Path, output_fasta: Path):
    """
    Runs MUSCLE to perform multiple sequence alignment.
    MUSCLE must be installed and accessible in your system's PATH.
    """
    try:
        # Ensure the output directory exists
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        # Run MUSCLE with input and output parameters
        subprocess.run(
            ["muscle", "-in", str(input_fasta), "-out", str(output_fasta)],
            check=True,
            stderr=subprocess.PIPE
        )
        print(f"Alignment completed. Aligned sequences saved in {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)

def calculate_conservation(alignment: MultipleSeqAlignment, threshold: float):
    """
    Calculates conservation for each column in the alignment.
    Returns two lists:
    - conserved_flags: List of boolean values indicating conserved positions.
    - most_common_aa: List of the most frequent amino acids in each position.
    """
    conserved_flags = []
    most_common_aa = []
    num_sequences = len(alignment)

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        # Exclude gaps and dots from conservation analysis
        column_no_gaps = column.replace('-', '').replace('.', '')
        if not column_no_gaps:
            conservation = 0
            most_common = '-'
        else:
            most_common = max(set(column_no_gaps), key=column_no_gaps.count)
            count = column_no_gaps.count(most_common)
            conservation = count / num_sequences
        conserved_flags.append(conservation >= threshold)
        most_common_aa.append(most_common)

    return conserved_flags, most_common_aa

def create_masked_alignments(alignment: MultipleSeqAlignment, conserved_flags):
    """
    Creates two masked alignments:
    1. Conserved: Keeps amino acids in conserved positions, replaces non-conserved with '-'.
    2. Non-Conserved: Keeps amino acids in non-conserved positions, replaces conserved with '-'.
    Returns two lists of SeqRecord.
    """
    conserved_alignment = []
    non_conserved_alignment = []

    for record in alignment:
        seq_conserved = []
        seq_non_conserved = []
        for i, is_conserved in enumerate(conserved_flags):
            residue = record.seq[i]
            if is_conserved:
                seq_conserved.append(residue)
                seq_non_conserved.append('-')
            else:
                seq_conserved.append('-')
                seq_non_conserved.append(residue)
        conserved_record = SeqRecord(Seq(''.join(seq_conserved)), id=record.id, description="")
        non_conserved_record = SeqRecord(Seq(''.join(seq_non_conserved)), id=record.id, description="")
        conserved_alignment.append(conserved_record)
        non_conserved_alignment.append(non_conserved_record)

    return conserved_alignment, non_conserved_alignment

def create_amalgam_sequence(alignment: MultipleSeqAlignment, conserved_flags, most_common_aa):
    """
    Creates an amalgam sequence:
    - If the position is conserved, inserts the most common amino acid.
    - Otherwise, inserts '-'.
    Returns a SeqRecord object.
    """
    amalgam_seq = []
    for i, is_conserved in enumerate(conserved_flags):
        if is_conserved:
            amalgam_seq.append(most_common_aa[i])
        else:
            amalgam_seq.append('-')
    amalgam_record = SeqRecord(
        Seq(''.join(amalgam_seq)),
        id="Amalgam_Sequence",
        description="Amalgamated sequence of aligned sequences"
    )
    return amalgam_record

def save_alignment(sequences, output_path):
    """
    Saves a list of SeqRecord to a FASTA file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(sequences, output_path, "fasta")
    print(f"Saved {len(sequences)} sequences in {output_path}")

def archive_sequences(sequences_path: Path, archive_path: Path):
    """
    Archives the sequences_path file into a ZIP archive.
    """
    try:
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(sequences_path, arcname=sequences_path.name)
        print(f"Archived {sequences_path} in {archive_path}")
    except Exception as e:
        print(f"Error archiving {sequences_path}: {e}", file=sys.stderr)

def main():
    # Step 1: Collect all sequences
    collect_sequences(DATA_DIR, COLLECTED_FASTA)

    # Step 2: Run MUSCLE for alignment
    run_muscle(COLLECTED_FASTA, ALIGNED_FASTA)

    # Step 3: Analyze the alignment
    try:
        alignment = AlignIO.read(ALIGNED_FASTA, "fasta")
        print(f"Alignment analyzed with {len(alignment)} sequences and {alignment.get_alignment_length()} positions.")
    except Exception as e:
        print(f"Error analyzing alignment: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 4: Calculate conservation
    conserved_flags, most_common_aa = calculate_conservation(alignment, CONSERVATION_THRESHOLD)
    num_conserved = sum(conserved_flags)
    print(f"Identified {num_conserved} conserved positions out of {len(conserved_flags)}.")

    # Step 5: Create masked alignments
    conserved_sequences, non_conserved_sequences = create_masked_alignments(alignment, conserved_flags)
    print(f"Masked alignments created: {len(conserved_sequences)} conserved sequences and {len(non_conserved_sequences)} non-conserved sequences.")

    # Step 6: Create the amalgam sequence
    amalgam_sequence = create_amalgam_sequence(alignment, conserved_flags, most_common_aa)
    print("Amalgam sequence created.")

    # Step 7: Save masked alignments and amalgam sequence
    conserved_fasta = RESULTS_DIR / "conserved_sequences.faa"
    non_conserved_fasta = RESULTS_DIR / "non_conserved_sequences.faa"
    amalgam_fasta = AMALGAM_FASTA

    save_alignment(conserved_sequences, conserved_fasta)
    save_alignment(non_conserved_sequences, non_conserved_fasta)
    save_alignment([amalgam_sequence], amalgam_fasta)

    # Step 8: Archive masked alignments and amalgam sequence
    archive_sequences(conserved_fasta, CONSERVE_ARCHIVE)
    archive_sequences(non_conserved_fasta, NON_CONSERVE_ARCHIVE)
    archive_sequences(amalgam_fasta, AMALGAM_ARCHIVE)

    print("All tasks completed successfully.")

if __name__ == "__main__":
    main()
