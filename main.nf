/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pacbio-hifi-variant-nextflow-pipeline — Multi-Sample Production Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Entry point for the PacBio HiFi variant calling pipeline.

    This file supports both single-sample and multi-sample (samplesheet) input.

    Responsibilities:
        1. Parameter validation and fail-fast checks
        2. Channel creation from samplesheet or single-sample input
        3. Calling the main workflow

    All pipeline logic lives in workflows/pacbio_variant_workflow.nf.

    Usage — single sample:
        nextflow run main.nf \
            --reads     test_data/synthetic_reads.fastq \
            --reference reference/chr20.fa \
            --sample    synthetic_reads \
            --outdir    results

    Usage — multi-sample samplesheet:
        nextflow run main.nf \
            --samplesheet samplesheet.csv \
            --reference   reference/chr20.fa \
            --outdir      results

    Samplesheet format (CSV, no header required):
        sample_name,/path/to/reads.fastq
        sample_A,data/sample_A.fastq
        sample_B,data/sample_B.fastq

    Usage — Docker:
        nextflow run main.nf -profile docker \
            --samplesheet samplesheet.csv \
            --reference   reference/chr20.fa \
            --outdir      results

    Test profile:
        nextflow run main.nf -profile test
        nextflow run main.nf -profile docker,test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { PACBIO_VARIANT_WORKFLOW } from './workflows/pacbio_variant_workflow'


// -----------------------------------------------------------------------
// Parameter defaults
// -----------------------------------------------------------------------
params.reads             = null
params.samplesheet       = null
params.reference         = null
params.sample            = null
params.outdir            = "results"
params.aligner           = "minimap2"
params.snv_caller        = "bcftools"
params.threads           = 4
params.clair3_model      = null
params.pipeline_container    = "pacbio-hifi-pipeline:1.0"
params.pbmm2_container       = "pacbio-hifi-pipeline:pbmm2"
params.clair3_container      = "pacbio-hifi-clair3:1.0"
params.deepvariant_container = "google/deepvariant:1.6.1"


// -----------------------------------------------------------------------
// Pipeline info
// -----------------------------------------------------------------------
log.info """
============================================
PacBio HiFi Variant Pipeline
============================================
samplesheet  : ${params.samplesheet ?: 'not provided'}
reads        : ${params.reads ?: 'not provided'}
reference    : ${params.reference}
outdir       : ${params.outdir}
aligner      : ${params.aligner}
snv_caller   : ${params.snv_caller}
threads      : ${params.threads}
============================================
""".stripIndent()


workflow {

    // -----------------------------------------------------------------------
    // SECTION 1 — Parameter validation
    // -----------------------------------------------------------------------

    // Reference is always required
    if (!params.reference) {
        error "Missing required parameter: --reference"
    }
    if (!file(params.reference).exists()) {
        error "Reference file not found: ${params.reference}"
    }

    // Either --samplesheet or --reads must be provided, not both
    if (!params.samplesheet && !params.reads) {
        error """
            No input provided. Use one of:
              --samplesheet samplesheet.csv   (multi-sample)
              --reads reads.fastq             (single sample, also requires --sample)
        """.stripIndent()
    }
    if (params.samplesheet && params.reads) {
        error "Provide either --samplesheet or --reads, not both."
    }

    // Single sample mode requires --sample
    if (params.reads && !params.sample) {
        error "Single sample mode requires --sample. Example: --sample my_sample"
    }

    // Single sample reads existence check
    if (params.reads && !file(params.reads).exists()) {
        error "Reads file not found: ${params.reads}"
    }

    // Samplesheet existence check
    if (params.samplesheet && !file(params.samplesheet).exists()) {
        error "Samplesheet not found: ${params.samplesheet}"
    }

    // Aligner validation
    def valid_aligners = ["minimap2", "pbmm2"]
    if (!valid_aligners.contains(params.aligner)) {
        error "Invalid aligner '${params.aligner}'. Valid options: ${valid_aligners.join(', ')}"
    }

    // Caller validation
    def valid_callers = ["bcftools", "clair3", "deepvariant"]
    if (!valid_callers.contains(params.snv_caller)) {
        error "Invalid snv_caller '${params.snv_caller}'. Valid options: ${valid_callers.join(', ')}"
    }

    // -----------------------------------------------------------------------
    // Docker-only tool checks
    // workflow.containerEngine is used instead of workflow.profile.contains()
    // because workflow.profile is not guaranteed to exist in all Nextflow
    // versions. containerEngine is set to 'docker', 'singularity', etc.
    // when a container engine is active, and null otherwise.
    // -----------------------------------------------------------------------
    if (params.aligner == "pbmm2" && !workflow.containerEngine) {
        error """
            pbmm2 requires -profile docker.
            Re-run with: nextflow run main.nf -profile docker --aligner pbmm2
            pbmm2 cannot be installed on Intel Mac via conda above v1.7.0.
        """.stripIndent()
    }
    if (params.snv_caller == "clair3" && !workflow.containerEngine) {
        error """
            clair3 requires -profile docker.
            Re-run with: nextflow run main.nf -profile docker --snv_caller clair3
        """.stripIndent()
    }
    if (params.snv_caller == "deepvariant" && !workflow.containerEngine) {
        error """
            deepvariant requires -profile docker.
            Re-run with: nextflow run main.nf -profile docker --snv_caller deepvariant
        """.stripIndent()
    }

    // Clair3 model validation
    if (params.snv_caller == "clair3" && !params.clair3_model) {
        error """
            clair3 requires a model directory: --clair3_model /path/to/model
            Download models from: https://github.com/HKU-BAL/Clair3#pre-trained-models
        """.stripIndent()
    }

    // Only check existence on host filesystem if path is not a container-internal path.
    // Container-internal paths (/opt/, /usr/) are valid when models are bundled in
    // the image. Host paths are validated here for early fail-fast behavior.
    if (params.snv_caller == "clair3" && params.clair3_model) {
        def model_path = params.clair3_model as String
        def container_paths = ["/opt/", "/usr/", "/home/"]
        def is_container_path = container_paths.any { model_path.startsWith(it) }
        if (!is_container_path && !file(params.clair3_model).exists()) {
            error "Clair3 model directory not found: ${params.clair3_model}"
        }
    }


    // -----------------------------------------------------------------------
    // SECTION 2 — Channel creation
    //
    // Supports two input modes:
    //
    // Mode A — Single sample (--reads + --sample)
    //   Creates a single-item channel from command line parameters.
    //
    // Mode B — Multi-sample samplesheet (--samplesheet)
    //   Reads a CSV file with format: sample_name,/path/to/reads.fastq
    //   Creates a channel with one item per sample.
    //   All samples flow through the pipeline in parallel.
    //
    // In both modes the channel structure is identical:
    //   [ val(sample), path(reads) ]
    // This ensures the workflow file requires no changes between modes.
    // -----------------------------------------------------------------------

    if (params.samplesheet) {

        // -----------------------------------------------------------------------
        // Multi-sample mode
        // splitCsv() parses each line of the CSV into a list.
        // map() converts each line into the standard tuple format.
        // Validates that each reads file exists before the pipeline starts.
        // -----------------------------------------------------------------------
        reads_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv()
            .map { row ->
                def sample = row[0].trim()
                def reads  = file(row[1].trim())

                if (!reads.exists()) {
                    error "Reads file not found for sample '${sample}': ${reads}"
                }
                if (sample.isEmpty()) {
                    error "Empty sample name found in samplesheet. Check your CSV."
                }

                tuple(sample, reads)
            }

        log.info "Multi-sample mode: reading from samplesheet ${params.samplesheet}"

    } else {

        // -----------------------------------------------------------------------
        // Single-sample mode
        // Channel.of() creates a single-item channel from command line parameters.
        // -----------------------------------------------------------------------
        reads_ch = Channel.of(
            tuple(params.sample, file(params.reads))
        )

        log.info "Single-sample mode: sample=${params.sample}"
    }

    ref_ch = Channel.fromPath(params.reference)


    // -----------------------------------------------------------------------
    // SECTION 3 — Workflow call
    // -----------------------------------------------------------------------

    PACBIO_VARIANT_WORKFLOW(reads_ch, ref_ch)
}


// -----------------------------------------------------------------------
// Completion handler
// -----------------------------------------------------------------------
workflow.onComplete {
    if (workflow.success) {
        log.info """
============================================
Pipeline completed successfully
============================================
Duration    : ${workflow.duration}
Output dir  : ${params.outdir}
Reports     : ${params.outdir}/reports/
============================================
""".stripIndent()
    } else {
        log.error """
============================================
Pipeline failed
============================================
Exit status : ${workflow.exitStatus}
Error msg   : ${workflow.errorMessage}
Work dir    : ${workflow.workDir}
============================================
""".stripIndent()
    }
}
