#!/usr/bin/env python3

"""
Generate synthetic PacBio HiFi reads from a reference region.

This script extracts a region from a reference FASTA and uses pbsim3
to simulate PacBio HiFi (CCS) reads for pipeline unit testing.

Usage example:
    python generate_synthetic_reads.py \
        --reference reference/chr20.fa \
        --chrom chr20 \
        --start 100000 \
        --end 118000 \
        --depth 20 \
        --outdir test_data

Note:
    The pbsim3 QSHMM-PACBIO-CCS model file must be downloaded separately.
    See README.md Setup Notes for instructions.

    A pre-generated test dataset for chr20:100000-118000 is provided
    in test_data/ for users who do not wish to regenerate reads.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd):
    """Run a shell command and exit with a clear message if it fails."""
    print(f"\n[CMD] {cmd}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        sys.exit(f"\nError: Command failed:\n  {cmd}\nCheck output above for details.")


def extract_region(reference, chrom, start, end, output_fasta):
    """
    Extract a genomic region from a FASTA file using seqkit.

    Args:
        reference:    Path to input FASTA file
        chrom:        Chromosome name (must match FASTA header exactly)
        start:        Start coordinate (1-based)
        end:          End coordinate (1-based, inclusive)
        output_fasta: Path to write extracted region
    """
    region = f"{start}:{end}"
    cmd = (
        f"seqkit subseq "
        f"--chr {chrom} "
        f"--region {region} "
        f"{reference} > {output_fasta}"
    )
    run_cmd(cmd)

    # Verify extraction produced output
    if not Path(output_fasta).exists() or Path(output_fasta).stat().st_size == 0:
        sys.exit(
            f"\nError: Region extraction produced no output.\n"
            f"Check that --chrom '{chrom}' matches the FASTA header exactly.\n"
            f"Check that --start {start} and --end {end} are within the sequence."
        )


def simulate_reads(reference_region, out_prefix, depth):
    """
    Simulate PacBio HiFi reads using pbsim3.

    Uses the QSHMM method with the PACBIO-CCS model to produce
    HiFi-like reads at Q20+ accuracy and ~15kb mean length.

    Args:
        reference_region: Path to extracted region FASTA
        out_prefix:       Output file prefix for pbsim3
        depth:            Target sequencing depth
    """
    cmd = (
	f"pbsim "
	f"--strategy wgs "
        f"--genome {reference_region} "
        f"--method qshmm "
	f"--qshmm $CONDA_PREFIX/data/QSHMM-PACBIO-CCS.model "
        f"--depth {depth} "
        f"--accuracy-mean 0.999 "
        f"--length-mean 15000 "
        f"--length-sd 2000 "
        f"--prefix {out_prefix} "
    )
    run_cmd(cmd)

"""
    Combine pbsim3 output FASTQ files into a single file.

    pbsim3 produces one FASTQ per reference sequence. This function
    merges them into a single file for downstream use.

    Args:
        out_prefix:   pbsim3 output prefix used during simulation
        output_fastq: Path to write combined FASTQ
"""  

def combine_reads(out_prefix, output_fastq):
    prefix_path = Path(out_prefix)
    
    # pbsim3 outputs .fq.gz or .fastq depending on version
    reads = sorted(prefix_path.parent.glob(f"{prefix_path.name}_*.fq.gz"))
    compressed = True
    
    if not reads:
        reads = sorted(prefix_path.parent.glob(f"{prefix_path.name}_*.fastq"))
        compressed = False

    if not reads:
        sys.exit(
            "\nError: No FASTQ files found after pbsim3 run.\n"
            "pbsim3 may have failed silently. Check logs above."
        )

    print(f"Combining {len(reads)} FASTQ file(s) (compressed={compressed})...")

    with open(output_fastq, "w") as outfile:
        for r in reads:
            if compressed:
                import gzip
                with gzip.open(r, "rt") as infile:
                    outfile.write(infile.read())
            else:
                with open(r) as infile:
                    outfile.write(infile.read())

"""
    Remove pbsim3 intermediate MAF and REF files.

    Args:
        out_prefix: pbsim3 output prefix used during simulation
"""

def cleanup_pbsim_intermediates(out_prefix):
    prefix_path = Path(out_prefix)
    removed = []

    for pattern in ["*.maf", "*.ref", "*.maf.gz", "*.fq.gz"]:
        for f in prefix_path.parent.glob(f"{prefix_path.name}_{pattern}"):
            f.unlink()
            removed.append(f.name)

    if removed:
        print(f"Removed intermediate files: {', '.join(removed)}")

def verify_output(fastq_path):
    """
    Print read statistics for the final FASTQ using seqkit.

    Args:
        fastq_path: Path to output FASTQ file
    """
    cmd = f"seqkit stats -a {fastq_path}"
    run_cmd(cmd)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate synthetic PacBio HiFi reads for pipeline testing.\n\n"
            "Extracts a region from a reference FASTA and simulates HiFi reads\n"
            "using pbsim3. Output is a FASTQ file ready for alignment with pbmm2."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA file (whole genome or single chromosome)"
    )

    parser.add_argument(
        "--chrom",
        required=True,
        help=(
            "Chromosome or sequence name to extract. "
            "Must match the FASTA header exactly (e.g. 'chr20' or '20'). "
            "Check your FASTA headers with: seqkit seq --name reference.fa"
        )
    )

    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Start coordinate of region to extract (1-based)"
    )

    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="End coordinate of region to extract (1-based, inclusive)"
    )

    parser.add_argument(
        "--depth",
        type=int,
        default=20,
        help="Target sequencing depth for simulation (default: 20)"
    )

    parser.add_argument(
        "--outdir",
        default="test_data",
        help="Output directory for generated files (default: test_data)"
    )

    args = parser.parse_args()

    # Validate region
    if args.end <= args.start:
        sys.exit(
            f"\nError: --end ({args.end}) must be greater than --start ({args.start})."
        )

    region_size = args.end - args.start
    if region_size < 15000:
        print(
            f"\nWarning: Region size ({region_size}bp) is smaller than the mean "
            f"read length (15000bp).\n"
            f"This may produce very few reads. Consider expanding your region\n"
            f"or reducing --length-mean in simulate_reads()."
        )

    # Set up output paths
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    region_fasta  = outdir / "reference_region.fa"
    out_prefix    = outdir / "synthetic_reads"
    final_fastq   = outdir / "synthetic_reads.fastq"

    # Run pipeline
    print(f"\nExtracting region {args.chrom}:{args.start}-{args.end}...")
    extract_region(args.reference, args.chrom, args.start, args.end, region_fasta)

    print("\nSimulating HiFi reads with pbsim3...")
    simulate_reads(region_fasta, out_prefix, args.depth)

    print("\nCombining reads...")
    combine_reads(out_prefix, final_fastq)

    print("\nCleaning up intermediate files...")
    cleanup_pbsim_intermediates(out_prefix)

    print("\nVerifying output...")
    verify_output(final_fastq)

    print("\n--- Done ---")
    print(f"Reference region : {region_fasta}")
    print(f"Synthetic reads  : {final_fastq}")
    print(
        f"\nExpected read count: ~{int(region_size * args.depth / 15000)} reads "
        f"({region_size}bp region at {args.depth}x depth, 15kb mean length)"
    )


if __name__ == "__main__":
    main()
