/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PACBIO_VARIANT_WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow for the PacBio HiFi variant calling pipeline.

    This file wires all modules together into a complete DAG (directed acyclic
    graph). It does not handle parameter parsing or channel creation — those
    are handled in main.nf. This workflow receives ready-made channels as input.

    Workflow:
        1. INDEX_REFERENCE     — index the reference FASTA
        2. ALIGN_*             — align HiFi reads (minimap2 or pbmm2)
        3. CALL_*              — call variants (bcftools, clair3, deepvariant)
        4. QC_SUMMARY          — generate QC report and HTML dashboard

    Channel wiring pattern:
        combine() is used to pair reads/BAM with the reference BEFORE branching.
        This is the correct nf-core pattern — it ensures every process receives
        a correctly structured single tuple and prevents channel race conditions
        in multi-sample runs.

    Aligner selection:
        --aligner minimap2     default, local + docker
        --aligner pbmm2        docker only

    Caller selection:
        --snv_caller bcftools      default, local + docker
        --snv_caller clair3        docker only
        --snv_caller deepvariant   docker only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// -----------------------------------------------------------------------
// Module imports
// Each module is imported from its file by process name.
// The 'include' statement makes the process available in this workflow.
// Nextflow resolves the file path relative to the pipeline root (projectDir).
//
// Local modules  — run with conda or docker
// Docker modules — require -profile docker
// -----------------------------------------------------------------------

include { INDEX_REFERENCE  } from '../modules/local/index_reference'
include { ALIGN_MINIMAP2   } from '../modules/local/align_minimap2'
include { CALL_BCFTOOLS    } from '../modules/local/call_bcftools'
include { QC_SUMMARY       } from '../modules/local/qc_summary'

include { ALIGN_PBMM2      } from '../modules/docker/align_pbmm2'
include { CALL_CLAIR3      } from '../modules/docker/call_clair3'
include { CALL_DEEPVARIANT } from '../modules/docker/call_deepvariant'


// -----------------------------------------------------------------------
// Workflow definition
//
// 'take' block — declares what this workflow receives as input channels.
// Channels are created in main.nf and passed in here.
// This separation means the workflow can be reused as a subworkflow
// in a larger pipeline without modification.
// -----------------------------------------------------------------------

workflow PACBIO_VARIANT_WORKFLOW {

    take:
    // reads_ch — tuple channel: [ val(sample), path(reads) ]
    // Each item represents one sample. For single-sample runs this
    // channel has one item. For multi-sample runs each sample flows
    // through all processes independently and in parallel.
    reads_ch

    // ref_ch — path channel: path(reference)
    // Single reference FASTA shared across all samples.
    ref_ch

    main:

    // -----------------------------------------------------------------------
    // STEP 1 — Index reference
    //
    // INDEX_REFERENCE takes the reference FASTA and produces a tuple:
    //   [ path(reference), path(fai) ]
    //
    // .out is how you access a process's output channel in Nextflow.
    //
    // .broadcast() ensures the reference channel can be consumed multiple
    // times by downstream processes without being exhausted after first use.
    // Without broadcast, the channel is consumed by the first process that
    // reads it and all subsequent processes block waiting for input that
    // never arrives — a common Nextflow pitfall in multi-process pipelines.
    // -----------------------------------------------------------------------

    INDEX_REFERENCE(ref_ch)
    ref_indexed_ch = INDEX_REFERENCE.out.first()


    // -----------------------------------------------------------------------
    // STEP 2 — Alignment
    //
    // combine() pairs every item in reads_ch with the indexed reference
    // BEFORE branching. This is the correct nf-core channel wiring pattern.
    //
    // Without combine(), passing two separate channels to a process causes
    // Nextflow to consume them independently — items may be mismatched in
    // multi-sample runs and the process input tuple will not match the
    // module declaration.
    //
    // After combine() each item in reads_ref has structure:
    //   [ val(sample), path(reads), path(reference), path(fai) ]
    //
    // branch() then routes each item to the correct aligner sub-channel.
    // Only the matching branch receives items — the other gets an empty
    // channel and Nextflow automatically skips processes with empty input.
    //
    // 'other' branch intentionally omitted — main.nf validates params.aligner
    // before the workflow runs so an invalid value can never reach this point.
    //
    // mix() reunifies BAM outputs from both aligners into a single channel
    // so downstream processes don't need to know which aligner was used.
    // -----------------------------------------------------------------------

    reads_ref = reads_ch.combine(ref_indexed_ch)

    reads_ref_branched = reads_ref.branch {
        minimap2: params.aligner == "minimap2"
        pbmm2:    params.aligner == "pbmm2"
    }

    ALIGN_MINIMAP2(reads_ref_branched.minimap2)
    ALIGN_PBMM2   (reads_ref_branched.pbmm2)

    // Reunify BAM outputs into a single channel regardless of aligner used
    bam_ch = ALIGN_MINIMAP2.out.mix(ALIGN_PBMM2.out)


    // -----------------------------------------------------------------------
    // STEP 3 — Variant calling
    //
    // Same combine-then-branch pattern as alignment.
    // bam_ref combines each BAM tuple with the indexed reference so every
    // caller receives a single structured channel item.
    //
    // After combine() each item in bam_ref has structure:
    //   [ val(sample), path(bam), path(bai), path(reference), path(fai) ]
    //
    // All three callers emit the same VCF tuple structure:
    //   [ val(sample), path(vcf), path(vcf_index) ]
    //
    // This consistent interface is what makes mix() work cleanly —
    // QC_SUMMARY doesn't need to know which caller produced the VCF.
    // -----------------------------------------------------------------------

    bam_ref = bam_ch.combine(ref_indexed_ch)

    bam_ref_branched = bam_ref.branch {
        bcftools:    params.snv_caller == "bcftools"
        clair3:      params.snv_caller == "clair3"
        deepvariant: params.snv_caller == "deepvariant"
    }

    CALL_BCFTOOLS   (bam_ref_branched.bcftools)
    CALL_CLAIR3     (bam_ref_branched.clair3)
    CALL_DEEPVARIANT(bam_ref_branched.deepvariant)

    // Reunify VCF outputs into a single channel regardless of caller used
    vcf_ch = CALL_BCFTOOLS.out
                .mix(CALL_CLAIR3.out)
                .mix(CALL_DEEPVARIANT.out)


    // -----------------------------------------------------------------------
    // STEP 4 — QC report
    //
    // QC_SUMMARY needs three inputs:
    //   1. BAM tuple    — [ val(sample), path(bam), path(bai) ]
    //   2. VCF tuple    — [ path(vcf), path(vcf_index) ]
    //   3. Reference    — [ path(reference), path(fai) ]
    //
    // join() matches items from bam_ch and vcf_ch by their first element
    // (index 0 = val(sample)). This is critical for multi-sample runs —
    // without joining, Nextflow could accidentally pair sample_A's BAM
    // with sample_B's VCF.
    //
    // After join() each item has structure:
    //   [ val(sample), path(bam), path(bai), path(vcf), path(vcf_idx) ]
    //
    // map() restructures this into the two separate tuples QC_SUMMARY
    // expects. The reference is passed separately as a shared channel.
    // -----------------------------------------------------------------------

    bam_vcf_ch = bam_ch
        .join(vcf_ch, by: 0)  // join on index 0 = sample name
        .map { sample, bam, bai, vcf, vcf_idx ->
            tuple(
                tuple(sample, bam, bai),  // BAM tuple for QC_SUMMARY input 1
                tuple(vcf, vcf_idx)       // VCF tuple for QC_SUMMARY input 2
            )
        }

    QC_SUMMARY(
        bam_vcf_ch.map { bam_tuple, vcf_tuple -> bam_tuple },
        bam_vcf_ch.map { bam_tuple, vcf_tuple -> vcf_tuple },
        ref_indexed_ch
    )

    // -----------------------------------------------------------------------
    // emit
    // Makes output channels available to main.nf or a parent workflow.
    // Not strictly required for a standalone pipeline but good practice —
    // it makes the workflow reusable as a subworkflow in a larger pipeline.
    // -----------------------------------------------------------------------

    emit:
    bam     = bam_ch
    vcf     = vcf_ch
    reports = QC_SUMMARY.out
}
