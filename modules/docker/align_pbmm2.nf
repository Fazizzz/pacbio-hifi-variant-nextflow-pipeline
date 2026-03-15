/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN_PBMM2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Aligns PacBio HiFi reads to a reference genome using pbmm2.

    pbmm2 is PacBio's native minimap2 wrapper optimized for PacBio read formats
    and chemistry presets. It produces better results than minimap2 on PacBio
    data due to tighter integration with PacBio read group metadata and
    chemistry-specific alignment parameters.

    DOCKER ONLY — pbmm2 cannot be installed on Intel Mac via conda above v1.7.0.
    This module requires -profile docker to run. The fail-fast check in main.nf
    catches any attempt to run pbmm2 without Docker before this process executes.

    Requires the pipeline container built with:
        docker build --build-arg ALIGNER=pbmm2 -t pacbio-hifi-pipeline:pbmm2 .

    Input:
        Combined tuple from workflow (reads_ch.combine(ref_indexed_ch)):
        [ val(sample), path(reads), path(reference), path(fai) ]

        fai is not passed as a flag — staged automatically by Nextflow
        into the work directory so pbmm2 can find it alongside the reference.

    Output:
        tuple (sample, bam, bai) — sorted, indexed BAM passed downstream

    Resource label:
        process_medium — multi-threaded alignment, moderate memory usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ALIGN_PBMM2 {

    // -----------------------------------------------------------------------
    // tag / label
    // -----------------------------------------------------------------------
    tag "${sample}"
    label 'process_medium'

    // -----------------------------------------------------------------------
    // conda — intentionally omitted
    // pbmm2 is not installable on Intel Mac via conda above v1.7.0.
    // This module is Docker only. The fail-fast check in main.nf catches
    // any attempt to run without -profile docker.
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // container
    // pbmm2 variant of the main pipeline image.
    // Built with: docker build --build-arg ALIGNER=pbmm2
    // Image version controlled via params.pbmm2_container in nextflow.config.
    // -----------------------------------------------------------------------
    container "${params.pbmm2_container}"

    // -----------------------------------------------------------------------
    // publishDir
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.{bam,bai}"

    input:
    // -----------------------------------------------------------------------
    // Combined tuple
    // Produced by reads_ch.combine(ref_indexed_ch) in the workflow.
    // path(fai) staged automatically — not passed as a flag argument.
    // pbmm2 locates reference.fa.fai by convention in the work directory.
    // -----------------------------------------------------------------------
    tuple val(sample), path(reads), path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    // -----------------------------------------------------------------------
    // pbmm2 vs minimap2 key differences:
    //
    // --preset HIFI
    //   pbmm2 uses named presets for chemistry-specific parameters.
    //   HIFI is equivalent to minimap2's map-hifi but with additional
    //   PacBio-specific optimizations for CCS read error profiles.
    //
    // --sort
    //   pbmm2 handles sorting natively in a single pass — no need to pipe
    //   through samtools sort. More memory efficient for large datasets.
    //
    // --unmapped
    //   Includes unmapped reads in the output BAM. Standard practice for
    //   PacBio data — maintains read count consistency with input.
    //
    // --rg
    //   Read group tag. pbmm2 uses --rg instead of minimap2's -R flag.
    //   Format is identical but the flag name differs.
    //
    // -j controls threading in pbmm2 (vs -t in minimap2)
    //
    // samtools index still required — pbmm2 sorts but does not index.
    // -----------------------------------------------------------------------
    """
    set -euo pipefail

    pbmm2 align \\
        --preset HIFI \\
        --sort \\
        --unmapped \\
        --rg "@RG\\tID:${sample}\\tSM:${sample}\\tPL:PACBIO" \\
        -j ${task.cpus} \\
        ${reference} \\
        ${reads} \\
        ${sample}.bam

    samtools index ${sample}.bam
    """
}
