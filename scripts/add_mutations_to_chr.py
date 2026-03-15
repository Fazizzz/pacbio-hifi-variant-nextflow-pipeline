#!/usr/bin/env python3

"""
Inject synthetic SNVs into a FASTA sequence and produce a truth VCF.

Designed for pipeline testing:

reference -> mutate -> simulate reads -> align -> call variants -> compare to truth

Supports:
- whole genome FASTA
- single chromosome FASTA
- mutation restriction to a specific region
"""

import argparse
import random
from pathlib import Path

BASES = ["A", "C", "G", "T"]


def read_fasta(fasta_path):
    """Load FASTA and return first sequence only."""
    headers = []
    seq = []
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if current_seq:
                    seq.append("".join(current_seq))
                    current_seq = []

                headers.append(line)

            else:
                current_seq.append(line)

    if current_seq:
        seq.append("".join(current_seq))

    if len(headers) > 1:
        print(f"Warning: {len(headers)} sequences found. Using first sequence only.")

    return headers[0], list(seq[0])


def write_fasta(header, seq, out_path, width=60):
    """Write FASTA with wrapped lines."""
    with open(out_path, "w") as f:

        f.write(header + "\n")

        for i in range(0, len(seq), width):
            f.write("".join(seq[i:i + width]) + "\n")


def choose_alt(ref_base):
    """Pick alternate base different from reference."""
    return random.choice([b for b in BASES if b != ref_base])


def inject_mutations(seq, chrom, n_mutations, seed, min_distance=100,
                     region_start=None, region_end=None):
    """
    Insert SNVs into sequence within optional region.
    Coordinates are 1-based genomic coordinates.
    """

    random.seed(seed)

    length = len(seq)

    if region_start is None:
        region_start = 1

    if region_end is None:
        region_end = length

    # Convert to Python 0-based indexing
    region_start = max(region_start - 1, 0)
    region_end = min(region_end - 1, length - 1)

    positions = []
    attempts = 0

    while len(positions) < n_mutations and attempts < 10000:

        pos = random.randint(region_start, region_end)

        if all(abs(pos - p) >= min_distance for p in positions):

            if seq[pos] in BASES:
                positions.append(pos)

        attempts += 1

    if len(positions) < n_mutations:
        print(
            f"Warning: Only able to place {len(positions)} mutations "
            f"with min_distance={min_distance}"
        )

    variants = []

    for pos in positions:

        ref = seq[pos]
        alt = choose_alt(ref)

        seq[pos] = alt

        variants.append({
            "chrom": chrom,
            "pos": pos + 1,
            "ref": ref,
            "alt": alt
        })

    return seq, variants


def write_vcf(variants, out_path, chrom, seq_length):
    """Write VCF containing injected variants."""

    # Ensure VCF is position sorted
    variants.sort(key=lambda v: v["pos"])

    with open(out_path, "w") as f:

        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={chrom},length={seq_length}>\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write(
            '##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,'
            'Description="Synthetically injected variant">\n'
        )

        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for v in variants:

            f.write(
                f"{v['chrom']}\t"
                f"{v['pos']}\t"
                f".\t"
                f"{v['ref']}\t"
                f"{v['alt']}\t"
                f".\t"
                f"PASS\tSYNTHETIC\n"
            )


def main():

    parser = argparse.ArgumentParser(
        description="Inject synthetic SNVs into a FASTA sequence."
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA file"
    )

    parser.add_argument(
        "--chrom",
        default="chr20",
        help="Chromosome name"
    )

    parser.add_argument(
        "--num_mutations",
        type=int,
        default=10,
        help="Number of SNVs to inject"
    )

    parser.add_argument(
        "--region_start",
        type=int,
        help="Start coordinate for mutation placement"
    )

    parser.add_argument(
        "--region_end",
        type=int,
        help="End coordinate for mutation placement"
    )

    parser.add_argument(
        "--min_distance",
        type=int,
        default=100,
        help="Minimum distance between mutations"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )

    parser.add_argument(
        "--outdir",
        default="test_data",
        help="Output directory"
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    print("Loading reference...")
    header, seq = read_fasta(args.reference)

    print("Injecting mutations...")

    seq, variants = inject_mutations(
        seq,
        args.chrom,
        args.num_mutations,
        args.seed,
        args.min_distance,
        args.region_start,
        args.region_end
    )

    mutated_fasta = outdir / "chr20_mutated.fa"
    truth_vcf = outdir / "truth_variants.vcf"

    print("Writing mutated FASTA...")
    write_fasta(header, seq, mutated_fasta)

    print("Writing truth VCF...")
    write_vcf(variants, truth_vcf, args.chrom, len(seq))

    print("\nMutation region:", args.region_start, "-", args.region_end)

    print("\nInjected", len(variants), "SNVs into", args.chrom)
    print("Positions:", [v["pos"] for v in variants])

    print("\nOutputs:")
    print("Mutated FASTA:", mutated_fasta)
    print("Truth VCF:", truth_vcf)


if __name__ == "__main__":
    main()
