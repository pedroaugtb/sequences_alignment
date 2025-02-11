#!/usr/bin/env python3
"""
This script is a modified version of your MAFFT alignment pipeline.
It now collects sequences only from genome folders in the 'complete' subdirectory,
runs MAFFT for multiple sequence alignment, calculates conservation, masks the alignment,
creates an amalgam sequence, and archives the results.
"""

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
# NOTE: DATA_DIR now points to the folder containing only complete genomes.
DATA_DIR = Path("/scratch/pedrobento/sequences_alignment/ncbi_dataset/ncbi_dataset/data/complete")
RESULTS_DIR = Path("results_mafft")
COLLECTED_FASTA = RESULTS_DIR / "all_sequences.faa"
ALIGNED_FASTA = RESULTS_DIR / "aligned.faa"
CONSERVE_ARCHIVE = RESULTS_DIR / "conserved_sequences.zip"
NON_CONSERVE_ARCHIVE = RESULTS_DIR / "non_conserved_sequences.zip"
AMALGAM_FASTA = RESULTS_DIR / "amalgam_sequence.faa"
AMALGAM_ARCHIVE = RESULTS_DIR / "amalgam_sequence.zip"
CONSERVATION_THRESHOLD = 0.2  # Increased threshold for stricter conservation

def collect_sequences(data_dir: Path, output_fasta: Path):
    """
    Collects all protein sequences from the genome folders (assumed to be complete)
    and writes them to a single FASTA file.
    """
    sequences = []
    missing_files = []  # To store folders where 'protein.faa' is missing

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

    if missing_files:
        print(f"Warning: {len(missing_files)} 'protein.faa' files were not found.", file=sys.stderr)
    else:
        print("All 'protein.faa' files were found.")

def run_mafft(input_fasta: Path, output_fasta: Path):
    """
    Runs MAFFT to perform a multiple sequence alignment.
    """
    try:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        with open(output_fasta, "w") as outfile:
            subprocess.run(
                ["mafft", "--auto", str(input_fasta)],
                check=True,
                stdout=outfile,
                stderr=subprocess.PIPE
            )
        print(f"Alignment completed. Aligned sequences saved in {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)

def calculate_conservation(alignment: MultipleSeqAlignment, threshold: float):
    """
    Calculates conservation for each column in the alignment.
    Returns two lists:
      - conserved_flags: Booleans indicating if the column is conserved.
      - most_common_aa: Most frequent amino acid in each column.
    """
    conserved_flags = []
    most_common_aa = []
    num_sequences = len(alignment)

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        # Remove gaps and dots
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
      1. Conserved: keeps amino acids in conserved positions (other positions are masked with '-').
      2. Non-Conserved: keeps amino acids in non-conserved positions (conserved positions masked).
    Returns two lists of SeqRecord objects.
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
    Creates an amalgam sequence where each position is filled with the most common amino acid if
    the position is conserved; otherwise, it is masked with '-'.
    Returns a SeqRecord.
    """
    amalgam_seq = []
    for i, is_conserved in enumerate(conserved_flags):
        amalgam_seq.append(most_common_aa[i] if is_conserved else '-')
    return SeqRecord(
        Seq(''.join(amalgam_seq)),
        id="Amalgam_Sequence",
        description="Amalgamated sequence of aligned sequences"
    )

def save_alignment(sequences, output_path):
    """
    Saves a list of SeqRecord objects to a FASTA file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(sequences, output_path, "fasta")
    print(f"Saved {len(sequences)} sequences in {output_path}")

def archive_sequences(sequences_path: Path, archive_path: Path):
    """
    Archives the given FASTA file into a ZIP archive.
    """
    try:
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(sequences_path, arcname=sequences_path.name)
        print(f"Archived {sequences_path} in {archive_path}")
    except Exception as e:
        print(f"Error archiving {sequences_path}: {e}", file=sys.stderr)

def main():
    # Step 1: Collect protein sequences from complete genome folders
    collect_sequences(DATA_DIR, COLLECTED_FASTA)

    # Step 2: Run MAFFT for multiple sequence alignment
    run_mafft(COLLECTED_FASTA, ALIGNED_FASTA)

    # Step 3: Read and analyze the alignment
    try:
        alignment = AlignIO.read(ALIGNED_FASTA, "fasta")
        print(f"Alignment loaded: {len(alignment)} sequences; {alignment.get_alignment_length()} positions.")
    except Exception as e:
        print(f"Error analyzing alignment: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 4: Calculate conservation metrics
    conserved_flags, most_common_aa = calculate_conservation(alignment, CONSERVATION_THRESHOLD)
    print(f"Identified {sum(conserved_flags)} conserved positions out of {len(conserved_flags)}.")

    # Step 5: Create masked alignments (conserved and non-conserved)
    conserved_sequences, non_conserved_sequences = create_masked_alignments(alignment, conserved_flags)
    print("Created masked alignments for conserved and non-conserved positions.")

    # Step 6: Create the amalgam sequence
    amalgam_sequence = create_amalgam_sequence(alignment, conserved_flags, most_common_aa)
    print("Amalgam sequence created.")

    # Step 7: Save results to FASTA files
    conserved_fasta = RESULTS_DIR / "conserved_sequences.faa"
    non_conserved_fasta = RESULTS_DIR / "non_conserved_sequences.faa"
    amalgam_fasta = AMALGAM_FASTA  # e.g., results_mafft/amalgam_sequence.faa

    save_alignment(conserved_sequences, conserved_fasta)
    save_alignment(non_conserved_sequences, non_conserved_fasta)
    save_alignment([amalgam_sequence], amalgam_fasta)

    # Step 8: Archive the FASTA files into ZIP archives
    archive_sequences(conserved_fasta, CONSERVE_ARCHIVE)
    archive_sequences(non_conserved_fasta, NON_CONSERVE_ARCHIVE)
    archive_sequences(amalgam_fasta, AMALGAM_ARCHIVE)

    print("All tasks completed successfully.")

if __name__ == "__main__":
    main()
