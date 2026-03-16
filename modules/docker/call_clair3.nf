/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL_CLAIR3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls small variants from a sorted PacBio HiFi BAM using Clair3.

    Clair3 is a deep-learning variant caller optimized for long reads.
    It produces higher sensitivity for HiFi data than bcftools, particularly
    for low-frequency variants and difficult genomic regions.

    DOCKER ONLY — Clair3 requires TensorFlow 2.x which creates irresolvable
    dependency conflicts with the main pipeline environment. It runs in a
    separate container based on hkubal/clair3:v0.1-r12.

    CONTAINER — built from containers/clair3/Dockerfile which installs
    bcftools and samtools on top of the hkubal/clair3 base image.
    bcftools is required for the normalization step after Clair3 inference.
    The base hkubal/clair3 image does not include bcftools.

    MODEL FILES — Clair3 requires a pre-trained model directory mounted at
    runtime. Models are not baked into the container (~1-2GB each).
    Download from: https://github.com/HKU-BAL/Clair3#pre-trained-models
    Pass the model directory via: --clair3_model /path/to/model

    Input:
        Combined tuple from workflow (bam_ch.combine(ref_indexed_ch)):
        [ val(sample), path(bam), path(bai), path(reference), path(fai) ]

        bai and fai staged automatically — not passed as flag arguments.

    Output:
        tuple (sample, vcf, vcf_index) — normalized VCF passed to QC_SUMMARY

    Resource label:
        process_high — deep learning inference is CPU/GPU intensive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CALL_CLAIR3 {

    // -----------------------------------------------------------------------
    // tag / label
    // process_high — Clair3 inference is significantly more resource
    // intensive than bcftools pileup calling.
    // -----------------------------------------------------------------------
    tag "${sample}"
    label 'process_high'

    // -----------------------------------------------------------------------
    // conda — intentionally omitted
    // Clair3 TensorFlow dependencies conflict with the main pipeline env.
    // Docker only. The fail-fast check in main.nf catches any attempt to
    // run without -profile docker.
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // container
    // Dedicated Clair3 image built from containers/clair3/Dockerfile.
    // Base: hkubal/clair3:v0.1-r12 (pinned — not latest).
    // bcftools installed in this image for the normalization step.
    // Image version controlled via params.clair3_container in nextflow.config.
    // -----------------------------------------------------------------------
    container "${params.clair3_container}"

    // -----------------------------------------------------------------------
    // publishDir
    // Only final normalized VCF and index published.
    // Clair3 produces many intermediate files — all stay in work/ only.
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: "*.vcf.gz*"

    input:
    // -----------------------------------------------------------------------
    // Combined tuple
    // Produced by bam_ch.combine(ref_indexed_ch) in the workflow.
    // bai and fai staged automatically — not passed as flag arguments.
    // -----------------------------------------------------------------------
    tuple val(sample), path(bam), path(bai), path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}_clair3_norm.vcf.gz"), path("${sample}_clair3_norm.vcf.gz.csi")

    script:
    // -----------------------------------------------------------------------
    // Model path validation
    // Checked before Clair3 inference so the pipeline fails fast with a
    // clear message rather than crashing with a Python stacktrace after
    // minutes of container startup and setup.
    //
    // Output VCF glob
    // Clair3 output naming varies across versions:
    //   merge_output.vcf.gz
    //   phased_merge_output.vcf.gz
    // Using a glob captures all variants safely.
    //
    // bcftools is installed in the Clair3 container via
    // containers/clair3/Dockerfile — not present in the base hkubal image.
    // -----------------------------------------------------------------------
    """
    set -euo pipefail

    # --------------------------------------------------
    # Validate model path before expensive inference
    # --------------------------------------------------

    if [ ! -d "${params.clair3_model}" ]; then
        echo "ERROR: Clair3 model directory not found: ${params.clair3_model}"
        echo "Download models from: https://github.com/HKU-BAL/Clair3#pre-trained-models"
        echo "Pass the model path with: --clair3_model /path/to/model"
        exit 1
    fi

    # --------------------------------------------------
    # Clair3 variant calling
    # --------------------------------------------------

    CLAIR3_OUTDIR=clair3_output

    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform="hifi" \\
        --model_path=${params.clair3_model} \\
        --output=\$CLAIR3_OUTDIR \\
        --include_all_ctgs \\
        --no_phasing_for_fa

    # --------------------------------------------------
    # Locate output VCF
    # Clair3 output naming varies across versions — glob captures all:
    #   merge_output.vcf.gz
    #   phased_merge_output.vcf.gz
    # --------------------------------------------------

    CLAIR3_VCF=\$(ls \$CLAIR3_OUTDIR/*merge_output*.vcf.gz 2>/dev/null | head -1)

    if [ -z "\$CLAIR3_VCF" ]; then
        echo "ERROR: No Clair3 output VCF found in \$CLAIR3_OUTDIR"
        echo "Contents of output directory:"
        ls -lh \$CLAIR3_OUTDIR/ || true
        exit 1
    fi

    echo "Using Clair3 output VCF: \$CLAIR3_VCF"

    # --------------------------------------------------
    # Normalize and rename to match pipeline conventions
    # bcftools installed in containers/clair3/Dockerfile
    # --------------------------------------------------

    bcftools norm \\
        --fasta-ref ${reference} \\
        --multiallelics -any \\
        --output-type z \\
        -o ${sample}_clair3_norm.vcf.gz \\
        \$CLAIR3_VCF

    bcftools index ${sample}_clair3_norm.vcf.gz
    """

    // -----------------------------------------------------------------------
    // stub
    // Used by -stub-run mode for lightweight CI validation.
    // Creates empty placeholder VCF and index so Nextflow can verify
    // channel wiring without executing Clair3 or Docker.
    // -----------------------------------------------------------------------
    stub:
    """
    touch ${sample}_clair3_norm.vcf.gz
    touch ${sample}_clair3_norm.vcf.gz.csi
    """
}
