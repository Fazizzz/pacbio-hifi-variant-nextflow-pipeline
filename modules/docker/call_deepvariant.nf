/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL_DEEPVARIANT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls small variants from a sorted PacBio HiFi BAM using DeepVariant.

    DeepVariant is a deep learning variant caller developed and maintained by
    Google. It is considered the highest-accuracy germline SNV/indel caller
    for PacBio HiFi data and is widely used in clinical and research settings
    at Illumina, PacBio, and most major sequencing centres.

    DOCKER ONLY — DeepVariant requires TensorFlow and specific system
    dependencies that conflict with the main pipeline environment.
    The fail-fast check in main.nf catches any attempt to run without Docker.

    NO CUSTOM DOCKERFILE — DeepVariant uses the official Google image directly.
    No build step required. The image is pulled automatically by Docker.
    Image pinned to a specific version — never use 'latest' in production
    pipelines as it breaks reproducibility.

    GPU SUPPORT — DeepVariant can use GPU acceleration if available.
    CPU mode is used by default. For GPU runs, add GPU flags to the
    docker.runOptions in nextflow.config.

    Input:
        Combined tuple from workflow (bam_ch.combine(ref_indexed_ch)):
        [ val(sample), path(bam), path(bai), path(reference), path(fai) ]

        bai and fai staged automatically — not passed as flag arguments.

    Output:
        tuple (sample, vcf, vcf_index) — normalized VCF passed to QC_SUMMARY
        gVCF is removed after normalization — not published downstream.

    Resource label:
        process_high — deep learning inference is CPU/GPU intensive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CALL_DEEPVARIANT {

    // -----------------------------------------------------------------------
    // tag / label
    // -----------------------------------------------------------------------
    tag "${sample}"
    label 'process_high'

    // -----------------------------------------------------------------------
    // conda — intentionally omitted
    // Docker only. The fail-fast check in main.nf catches any attempt to
    // run without -profile docker.
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // container
    // Official Google DeepVariant image — no custom build required.
    // Pinned to a specific version for reproducibility.
    // Image version controlled via params.deepvariant_container in
    // nextflow.config so it can be updated in one place.
    // -----------------------------------------------------------------------
    container "${params.deepvariant_container}"

    // -----------------------------------------------------------------------
    // publishDir
    // Only final normalized VCF and index are published.
    // gVCF is removed in the script block before publishDir is evaluated
    // so it is never copied to the output directory.
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: "*.vcf.gz*"

    input:
    // -----------------------------------------------------------------------
    // Combined tuple
    // Produced by bam_ch.combine(ref_indexed_ch) in the workflow.
    // bai and fai staged automatically — not passed as flag arguments.
    // DeepVariant locates .bai and .fai by convention in the work directory.
    // -----------------------------------------------------------------------
    tuple val(sample), path(bam), path(bai), path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}_deepvariant_norm.vcf.gz"), path("${sample}_deepvariant_norm.vcf.gz.csi")

    script:
    // -----------------------------------------------------------------------
    // /opt/deepvariant/bin/run_deepvariant
    //   Explicit binary path — not guaranteed to be on PATH across all
    //   DeepVariant image versions. Explicit path is consistent across
    //   all official Google DeepVariant images.
    //
    // --model_type PACBIO
    //   Uses PacBio HiFi-specific model weights bundled in the container.
    //   Do not use WGS or WES model types for HiFi data.
    //   No external model files required — unlike Clair3.
    //
    // --num_shards
    //   Controls parallelism in the make_examples stage.
    //   Set to task.cpus for maximum throughput.
    //
    // gVCF cleanup
    //   DeepVariant always produces a gVCF alongside the VCF.
    //   Deleted after normalization to keep results/ clean.
    //   The publishDir pattern *.vcf.gz* would otherwise capture it.
    //   To enable gVCF output in future, remove the rm line and add
    //   the gVCF path to the output block.
    //
    // bcftools norm
    //   Ensures consistent VCF format for comparison with bcftools
    //   and Clair3 outputs and downstream benchmarking with RTG vcfeval.
    //   bcftools is available in the official DeepVariant image.
    // -----------------------------------------------------------------------
    """
    set -euo pipefail

    RAW_VCF=${sample}_deepvariant_raw.vcf.gz
    GVCF=${sample}_deepvariant.g.vcf.gz
    FINAL_VCF=${sample}_deepvariant_norm.vcf.gz

    # --------------------------------------------------
    # DeepVariant variant calling
    # Three internal stages handled by run_deepvariant:
    #   1. make_examples  — pileup image generation
    #   2. call_variants  — neural network inference
    #   3. postprocess    — VCF output
    # --------------------------------------------------

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=PACBIO \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=\$RAW_VCF \\
        --output_gvcf=\$GVCF \\
        --num_shards=${task.cpus}

    # --------------------------------------------------
    # Remove gVCF — not published downstream
    # Delete before publishDir evaluates to prevent
    # the *.vcf.gz* glob from capturing it.
    # --------------------------------------------------

    rm -f \$GVCF \${GVCF}.tbi \${GVCF}.csi 2>/dev/null || true

    # --------------------------------------------------
    # Normalize to match pipeline VCF conventions
    # --------------------------------------------------

    bcftools norm \\
        --fasta-ref ${reference} \\
        --multiallelics -any \\
        --output-type z \\
        -o \$FINAL_VCF \\
        \$RAW_VCF

    bcftools index \$FINAL_VCF
    """
}
