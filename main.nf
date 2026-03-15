/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pacbio-hifi-variant-nextflow-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Entry point for the PacBio HiFi variant calling pipeline.

    This file is deliberately minimal. Its only responsibilities are:
        1. Parameter validation and fail-fast checks
        2. Channel creation from user inputs
        3. Calling the main workflow

    All pipeline logic lives in workflows/pacbio_variant_workflow.nf.

    Usage:
        # Default run — minimap2 + bcftools, local conda environment
        nextflow run main.nf \
            --reads     test_data/synthetic_reads.fastq \
            --reference reference/chr20.fa \
            --sample    synthetic_reads \
            --outdir    results

        # Docker run — minimap2 + bcftools in container
        nextflow run main.nf \
            -profile docker \
            --reads     test_data/synthetic_reads.fastq \
            --reference reference/chr20.fa \
            --sample    synthetic_reads \
            --outdir    results

        # Docker run — pbmm2 aligner
        nextflow run main.nf \
            -profile docker \
            --aligner   pbmm2 \
            --reads     test_data/synthetic_reads.fastq \
            --reference reference/chr20.fa \
            --sample    synthetic_reads

        # Docker run — Clair3 variant caller
        nextflow run main.nf \
            -profile    docker \
            --snv_caller clair3 \
            --clair3_model /path/to/model \
            --reads     test_data/synthetic_reads.fastq \
            --reference reference/chr20.fa \
            --sample    synthetic_reads

        # Test profile — runs on provided synthetic chr20 test dataset
        nextflow run main.nf -profile test
        nextflow run main.nf -profile docker,test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// -----------------------------------------------------------------------
// Nextflow DSL version
// DSL2 is the modern Nextflow syntax used throughout this pipeline.
// It enables module imports, named workflows, and channel operations
// like branch(), join(), mix(), combine(), and broadcast().
// -----------------------------------------------------------------------
nextflow.enable.dsl = 2

// -----------------------------------------------------------------------
// Workflow import
// Imports the named workflow from the workflow file.
// This is the only import in main.nf — all module imports are handled
// inside pacbio_variant_workflow.nf.
// -----------------------------------------------------------------------
include { PACBIO_VARIANT_WORKFLOW } from './workflows/pacbio_variant_workflow'


// -----------------------------------------------------------------------
// Parameter defaults
// All pipeline parameters declared here with their default values.
// Users override these on the command line with --param_name value.
// Parameters can also be set in params.yaml passed with -params-file.
//
// null means the parameter is required — the validation block below
// catches any missing required parameters before the pipeline runs.
// -----------------------------------------------------------------------
params.reads                 = null
params.reference             = null
params.sample                = null
params.outdir                = "results"
params.aligner               = "minimap2"
params.snv_caller            = "bcftools"
params.threads               = 4
params.clair3_model          = null
params.pipeline_container    = "pacbio-hifi-pipeline:1.0"
params.pbmm2_container       = "pacbio-hifi-pipeline:pbmm2"
params.clair3_container      = "pacbio-hifi-clair3:1.0"
params.deepvariant_container = "google/deepvariant:1.6.1"


// -----------------------------------------------------------------------
// Pipeline info block
// Prints a summary header to the log at pipeline startup.
// log.info writes to the Nextflow log — visible in terminal and
// in the .nextflow.log file for debugging.
// -----------------------------------------------------------------------
log.info """
============================================
PacBio HiFi Variant Pipeline
============================================
reads        : ${params.reads}
reference    : ${params.reference}
sample       : ${params.sample}
outdir       : ${params.outdir}
aligner      : ${params.aligner}
snv_caller   : ${params.snv_caller}
threads      : ${params.threads}
============================================
""".stripIndent()


// -----------------------------------------------------------------------
// Main workflow
// Three sections:
//   1. Parameter validation — fail fast with clear messages
//   2. Channel creation     — convert file paths to Nextflow channels
//   3. Workflow call        — pass channels to PACBIO_VARIANT_WORKFLOW
// -----------------------------------------------------------------------
workflow {

    // -----------------------------------------------------------------------
    // SECTION 1 — Parameter validation
    //
    // All validation happens here before any processes run.
    // 'error' immediately stops the pipeline with a clear message.
    // This is preferable to letting tools fail with cryptic errors
    // after minutes of setup.
    // -----------------------------------------------------------------------

    // Required parameters
    if (!params.reads) {
        error "Missing required parameter: --reads"
    }
    if (!params.reference) {
        error "Missing required parameter: --reference"
    }
    if (!params.sample) {
        error "Missing required parameter: --sample"
    }

    // Input file existence checks
    if (!file(params.reads).exists()) {
        error "Reads file not found: ${params.reads}"
    }
    if (!file(params.reference).exists()) {
        error "Reference file not found: ${params.reference}"
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
    //
    // workflow.containerEngine is used instead of workflow.profile.contains()
    // because workflow.profile is not guaranteed to exist in all Nextflow
    // versions. containerEngine is set to 'docker', 'singularity', etc.
    // when a container engine is active, and null otherwise — making it
    // a reliable cross-version check.
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

    // Clair3 model path validation
    // Checked here so the error is caught before any Docker container starts
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
    // Channels are the pipes that carry data between processes.
    // Here we convert user-supplied file paths into Nextflow channels.
    //
    // Channel.of() creates a channel from one or more values.
    //
    // tuple(params.sample, file(params.reads))
    //   Creates a tuple channel item pairing the sample name with the
    //   reads file. The sample name travels with its data throughout
    //   the entire pipeline so outputs are always named correctly.
    //   This is the standard nf-core pattern for sample-aware pipelines.
    //
    // Channel.fromPath()
    //   Creates a channel from a file path. Used for the reference since
    //   it is shared across all samples and doesn't carry a sample name.
    // -----------------------------------------------------------------------

    reads_ch = Channel.of(
        tuple(params.sample, file(params.reads))
    )

    ref_ch = Channel.fromPath(params.reference)


    // -----------------------------------------------------------------------
    // SECTION 3 — Workflow call
    //
    // Pass the channels to the main workflow.
    // Everything else — module imports, process execution, channel wiring —
    // is handled inside pacbio_variant_workflow.nf.
    // -----------------------------------------------------------------------

    PACBIO_VARIANT_WORKFLOW(reads_ch, ref_ch)
}


// -----------------------------------------------------------------------
// Completion handler
// Runs after the pipeline finishes — prints a summary to the log.
// workflow.success is a Nextflow built-in boolean.
// workflow.duration is automatically tracked by Nextflow.
// Runs regardless of success or failure so you always get a summary.
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
