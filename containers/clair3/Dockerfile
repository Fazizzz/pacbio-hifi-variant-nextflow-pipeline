FROM hkubal/clair3:v0.1-r12

LABEL org.opencontainers.image.title="PacBio HiFi Variant Pipeline - Clair3"
LABEL org.opencontainers.image.description="Clair3 variant calling container for the PacBio HiFi Nextflow pipeline"
LABEL org.opencontainers.image.version="1.0"

# -----------------------------------------------------------------------
# Base image: hkubal/clair3:v0.1-r12 (pinned — not latest)
# Includes: clair3, TensorFlow, run_clair3.sh entrypoint
# Base runs Python 3.9 on glibc 2.31 (Linux/amd64)
# Does NOT include: bcftools, samtools (added below)
#
# bcftools and samtools are required for the normalization step
# that runs after Clair3 inference in call_clair3.nf.
#
# Version pins removed for samtools and bcftools — the Clair3 base
# image runs Python 3.9 on glibc 2.31 which conflicts with
# samtools/bcftools 1.19 system library requirements.
# Conda resolves the highest compatible versions automatically.
# -----------------------------------------------------------------------

USER root

# -----------------------------------------------------------------------
# Install post-processing utilities
# bioconda before conda-forge for correct package resolution.
# No version pins — let conda resolve compatible versions with
# the Clair3 base environment.
# -----------------------------------------------------------------------
RUN conda install -y \
    -c bioconda \
    -c conda-forge \
    samtools \
    bcftools && \
    conda clean -a

WORKDIR /pipeline

# -----------------------------------------------------------------------
# Smoke test
# --help used instead of --version for run_clair3.sh —
# --version is not a guaranteed flag across Clair3 releases.
# Print installed versions for documentation purposes.
# -----------------------------------------------------------------------
RUN run_clair3.sh --help 2>&1 | head -5 && \
    samtools --version | head -1 && \
    bcftools --version | head -1

CMD ["bash"]
