/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN_MINIMAP2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Aligns PacBio HiFi reads to a reference genome using minimap2.

    This is the default aligner for this pipeline. It runs locally via conda
    and inside the main pipeline Docker container.

    Input:
        Combined tuple from workflow (reads_ch.combine(ref_indexed_ch)):
        [ val(sample), path(reads), path(reference), path(fai) ]

        The workflow combines reads_ch and ref_indexed_ch before calling
        this module so all inputs arrive as a single structured channel item.
        fai is not passed as a flag — it is staged automatically by Nextflow
        into the work directory so minimap2 can find it alongside the reference.

    Output:
        tuple (sample, bam, bai) — sorted, indexed BAM passed downstream.
        sample name travels with the BAM so downstream processes name
        outputs correctly in multi-sample runs.

    Resource label:
        process_medium — multi-threaded alignment, moderate memory usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ALIGN_MINIMAP2 {

    // -----------------------------------------------------------------------
    // tag
    // Uses the sample name so each alignment job is clearly identified
    // in logs and the Nextflow execution report.
    // "${sample}" syntax preferred over "$sample" in directives
    // to avoid binding issues when Nextflow resolves the value.
    // -----------------------------------------------------------------------
    tag "${sample}"

    // -----------------------------------------------------------------------
    // label
    // process_medium — defined in nextflow.config
    // Typically 4 CPUs, 8 GB memory for alignment workloads.
    // -----------------------------------------------------------------------
    label 'process_medium'

    // -----------------------------------------------------------------------
    // conda
    // Packages installed when running with -profile conda or locally.
    // Nextflow creates and caches the environment automatically.
    // bioconda listed first for correct package resolution.
    // -----------------------------------------------------------------------
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.19"

    // -----------------------------------------------------------------------
    // container
    // References a param defined in nextflow.config rather than hardcoding
    // the image name. Update image version once in nextflow.config without
    // touching module files.
    //
    // params.pipeline_container = "pacbio-hifi-pipeline:1.0"
    // -----------------------------------------------------------------------
    container "${params.pipeline_container}"

    // -----------------------------------------------------------------------
    // publishDir
    // Sorted BAM and index published to outdir/bam/
    // pattern: restricts publishing to only .bam and .bai files.
    // Intermediate files produced during the script stay in work/ only.
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.{bam,bai}"

    input:
    // -----------------------------------------------------------------------
    // Combined tuple
    // Produced by reads_ch.combine(ref_indexed_ch) in the workflow.
    // All four values arrive together as one channel item — this is the
    // correct nf-core pattern for pairing sample data with a shared reference.
    //
    // path(fai) — not referenced in the script as a flag argument.
    // Nextflow stages it into the work directory automatically so minimap2
    // can find reference.fa.fai alongside reference.fa by convention.
    // -----------------------------------------------------------------------
    tuple val(sample), path(reads), path(reference), path(fai)

    output:
    // -----------------------------------------------------------------------
    // tuple output
    // sample name + BAM + index travel together as one channel item.
    // Downstream processes (CALL_BCFTOOLS, QC_SUMMARY) destructure this
    // tuple to access each component by name.
    // -----------------------------------------------------------------------
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    // -----------------------------------------------------------------------
    // set -euo pipefail
    // Nextflow script blocks do not inherit strict mode automatically.
    // Without this, a failure in the minimap2 | samtools pipe could
    // go undetected and produce an empty or corrupt BAM silently.
    //
    // task.cpus
    // Injected by Nextflow at runtime from the process_medium label
    // defined in nextflow.config. Never hardcode thread counts in modules —
    // use task.cpus so resources can be tuned in one place.
    //
    // --secondary=no
    // Suppresses secondary alignments — only primary alignments are kept.
    // Secondary reads add noise to variant callers.
    //
    // -R read group tag
    // bcftools and most variant callers expect a read group in the BAM header.
    // -----------------------------------------------------------------------
    """
    set -euo pipefail

    minimap2 \\
        -t ${task.cpus} \\
        -ax map-hifi \\
        --secondary=no \\
        -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:PACBIO" \\
        ${reference} \\
        ${reads} \\
    | samtools view -b - \\
    | samtools sort -@ ${task.cpus} -o ${sample}.bam

    samtools index ${sample}.bam
    """

    // -----------------------------------------------------------------------
    // stub
    // Used by -stub-run mode for lightweight CI validation.
    // Creates empty placeholder BAM and index files so Nextflow can verify
    // channel wiring without executing minimap2 or samtools.
    // -----------------------------------------------------------------------
    stub:
    """
    touch ${sample}.bam
    touch ${sample}.bam.bai
    """
}
