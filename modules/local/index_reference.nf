/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEX_REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generates a samtools FASTA index (.fai) for the reference genome.

    This process runs once at pipeline startup and passes the indexed
    reference to all downstream processes that require it.

    Input:
        path reference  —  reference FASTA file

    Output:
        tuple of (reference, fai) passed as a single channel item
        so downstream processes always receive both files together.

    Resource label:
        process_low — single-threaded, lightweight operation.
        Labels are defined in nextflow.config and allow resource tuning
        in one place rather than editing every module individually.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INDEX_REFERENCE {

    // -----------------------------------------------------------------------
    // tag
    // A label shown in Nextflow logs for this process execution.
    // Makes it easy to identify which file is being processed
    // when reading logs or the Nextflow execution report.
    // -----------------------------------------------------------------------
    tag "$reference"

    // -----------------------------------------------------------------------
    // label
    // Links this process to a resource profile defined in nextflow.config.
    // 'process_low' = 1 CPU, 1 GB memory — appropriate for samtools faidx
    // which is single-threaded and completes in seconds.
    //
    // Other labels used in this pipeline:
    //   process_medium — alignment (minimap2, pbmm2)
    //   process_high   — variant calling (clair3, deepvariant)
    // -----------------------------------------------------------------------
    label 'process_low'

    // -----------------------------------------------------------------------
    // cpus / memory
    // Explicit resource declarations for this process.
    // These override the label defaults if both are present,
    // but here they match process_low for clarity and documentation.
    // -----------------------------------------------------------------------
    cpus   1
    memory '1 GB'

    // -----------------------------------------------------------------------
    // publishDir
    // Copies output files to a user-facing directory after the process
    // completes. Without this, outputs stay buried in work/ and are
    // inaccessible to the user.
    //
    // params.outdir is set by --outdir on the command line (default: results)
    // mode: 'copy' copies files rather than symlinking (safer for portability)
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/reference", mode: 'copy'

    // -----------------------------------------------------------------------
    // input
    // Declares what this process receives from a channel.
    // 'path' tells Nextflow this is a file — it will be automatically
    // staged (linked) into the process work directory before the
    // script runs. You do not need to handle file paths manually.
    // -----------------------------------------------------------------------
    input:
    path reference

    // -----------------------------------------------------------------------
    // output
    // Declares what this process emits back into a channel.
    //
    // 'tuple' groups multiple values into a single channel item.
    // This is the standard Nextflow pattern for keeping related files
    // together as they flow through the pipeline. Downstream processes
    // that need both the reference and its index receive them as one unit.
    //
    // '*.fai' glob captures the index regardless of the reference filename.
    // Each process runs in an isolated work directory so there is no
    // risk of capturing unrelated .fai files.
    // -----------------------------------------------------------------------
    output:
    tuple path(reference), path("*.fai")

    // -----------------------------------------------------------------------
    // script
    // The shell commands to run inside the process work directory.
    // Nextflow variables are interpolated with $ syntax.
    // Each process execution gets its own clean isolated directory
    // under work/ — this is how Nextflow achieves reproducibility
    // and enables the -resume flag to skip completed steps.
    // -----------------------------------------------------------------------
    script:
    """
    samtools faidx ${reference}
    """

    // -----------------------------------------------------------------------
    // stub
    // Used by -stub-run mode for lightweight CI validation.
    // Creates empty placeholder files matching the declared output types
    // so Nextflow can verify channel wiring without executing real tools.
    // -----------------------------------------------------------------------
    stub:
    """
    touch ${reference}.fai
    """
}
