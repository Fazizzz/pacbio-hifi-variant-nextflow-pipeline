/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL_BCFTOOLS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls small variants from a sorted PacBio HiFi BAM using bcftools.

    Default SNV caller for this pipeline. Runs locally via conda and inside
    the main pipeline Docker container — no additional setup required.

    Pipeline:
        mpileup → call → raw VCF
                       → filter (QUAL + INFO/DP, soft-filter LowQual)
                       → normalize
                       → final VCF

    Input:
        Combined tuple from workflow (bam_ch.combine(ref_indexed_ch)):
        [ val(sample), path(bam), path(bai), path(reference), path(fai) ]

        bai and fai are not passed as flag arguments — staged automatically
        by Nextflow so bcftools can find them alongside their parent files
        by convention in the work directory.

    Output:
        tuple (sample, vcf, vcf_index) — normalized VCF passed to QC_SUMMARY

    Resource label:
        process_medium — CPU-moderate, memory light
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CALL_BCFTOOLS {

    // -----------------------------------------------------------------------
    // tag / label
    // -----------------------------------------------------------------------
    tag "${sample}"
    label 'process_medium'

    // -----------------------------------------------------------------------
    // conda
    // Packages installed when running with -profile conda or locally.
    // -----------------------------------------------------------------------
    conda "bioconda::bcftools=1.19 bioconda::samtools=1.19"

    // -----------------------------------------------------------------------
    // container
    // Same main pipeline image as ALIGN_MINIMAP2 — bcftools and samtools
    // are included in pacbio-hifi-pipeline:1.0.
    // Image version controlled via params.pipeline_container in nextflow.config.
    // -----------------------------------------------------------------------
    container "${params.pipeline_container}"

    // -----------------------------------------------------------------------
    // publishDir
    // Only final VCF and index published to outdir/vcf/
    // pattern: "*.vcf.gz*" matches both .vcf.gz and .vcf.gz.csi
    // Intermediate raw and filtered VCFs stay in work/ only.
    // -----------------------------------------------------------------------
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: "*.vcf.gz*"

    input:
    // -----------------------------------------------------------------------
    // Combined tuple
    // Produced by bam_ch.combine(ref_indexed_ch) in the workflow.
    // bai staged automatically — bcftools requires it alongside the BAM.
    // fai staged automatically — bcftools mpileup requires it alongside
    // the reference. Neither is passed as an explicit flag argument.
    // -----------------------------------------------------------------------
    tuple val(sample), path(bam), path(bai), path(reference), path(fai)

    output:
    // -----------------------------------------------------------------------
    // sample name travels with the VCF so QC_SUMMARY and any downstream
    // processes know which sample produced these variants.
    // -----------------------------------------------------------------------
    tuple val(sample), path("${sample}_norm.vcf.gz"), path("${sample}_norm.vcf.gz.csi")

    script:
    // -----------------------------------------------------------------------
    // Variable escaping in Nextflow script blocks:
    //   ${sample}    — Nextflow variable, interpolated before script runs
    //   \$RAW_VCF    — bash variable, backslash escapes $ so Nextflow
    //                  does not interpret it as a Nextflow variable
    //
    // task.cpus — injected at runtime from the process_medium label.
    // Never hardcode thread counts in modules.
    //
    // INFO/DP vs FORMAT/DP — bcftools mpileup with FORMAT/AD,FORMAT/DP
    // annotation creates both INFO/DP and FORMAT/DP in the VCF header.
    // The filter must explicitly reference INFO/DP (site-level depth)
    // to avoid an ambiguous filtering expression error.
    //
    // --soft-filter LowQual — tags failing variants in the FILTER column
    // rather than silently dropping them. Makes the VCF auditable.
    //
    // BQSR not used — BQSR is a GATK short-read step designed to correct
    // systematic Illumina base quality errors. PacBio HiFi reads use a
    // different error model and are already quality-calibrated by the
    // CCS consensus process. Applying BQSR to HiFi data is incorrect.
    // -----------------------------------------------------------------------
    """
    set -euo pipefail

    RAW_VCF=${sample}_raw.vcf.gz
    FILTERED_VCF=${sample}_filtered.vcf.gz
    FINAL_VCF=${sample}_norm.vcf.gz

    # --------------------------------------------------
    # Variant calling
    # mpileup generates per-base pileup, bcftools call
    # performs genotyping using the multiallelic model.
    # --------------------------------------------------

    bcftools mpileup \\
        --fasta-ref ${reference} \\
        --output-type u \\
        --max-depth 200 \\
        --min-MQ 20 \\
        --min-BQ 20 \\
        --annotate FORMAT/AD,FORMAT/DP \\
        --threads ${task.cpus} \\
        ${bam} \\
    | bcftools call \\
        --multiallelic-caller \\
        --variants-only \\
        --output-type z \\
        --threads ${task.cpus} \\
        -o \$RAW_VCF

    bcftools index \$RAW_VCF

    # --------------------------------------------------
    # Filtering
    # INFO/DP used explicitly to avoid ambiguity with
    # FORMAT/DP which is also present in the VCF header.
    # LowQual soft-filter preserves filtered records
    # for audit — they are not silently dropped.
    # --------------------------------------------------

    bcftools filter \\
        --include 'QUAL >= 20 && INFO/DP >= 5' \\
        --soft-filter LowQual \\
        --output-type z \\
        -o \$FILTERED_VCF \\
        \$RAW_VCF

    bcftools index \$FILTERED_VCF

    # --------------------------------------------------
    # Normalization
    # Left-aligns and splits multiallelic sites.
    # Required for correct comparison against truth VCF
    # using downstream benchmarking tools (e.g. RTG vcfeval).
    # --------------------------------------------------

    bcftools norm \\
        --fasta-ref ${reference} \\
        --multiallelics -any \\
        --output-type z \\
        -o \$FINAL_VCF \\
        \$FILTERED_VCF

    bcftools index \$FINAL_VCF
    """

    // -----------------------------------------------------------------------
    // stub
    // Used by -stub-run mode for lightweight CI validation.
    // Creates empty placeholder VCF and index files so Nextflow can verify
    // channel wiring without executing bcftools.
    // -----------------------------------------------------------------------
    stub:
    """
    touch ${sample}_norm.vcf.gz
    touch ${sample}_norm.vcf.gz.csi
    """
}
