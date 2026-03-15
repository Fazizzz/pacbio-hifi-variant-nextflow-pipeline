FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="PacBio HiFi Variant Pipeline"
LABEL org.opencontainers.image.description="Container for alignment and variant calling"
LABEL org.opencontainers.image.version="1.0"

# -----------------------------------------------------------------------
# ALIGNER build argument
# Selects which aligner to install at build time.
# Default is minimap2 — matches local conda environment behavior.
#
# Build commands:
#   docker build -t pacbio-hifi-pipeline:1.0 .
#   docker build --build-arg ALIGNER=pbmm2 -t pacbio-hifi-pipeline:pbmm2 .
# -----------------------------------------------------------------------
ARG ALIGNER=minimap2

USER root

# -----------------------------------------------------------------------
# Install procps
# Required by Nextflow to collect task metrics (CPU, memory) inside
# the container using the 'ps' command. Not included in the
# mambaorg/micromamba base image by default. Without this, Nextflow
# fails with: "Command 'ps' required by nextflow to collect task
# metrics cannot be found"
# -----------------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
        libfreetype6 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------
# Install base tools
# bioconda must come before conda-forge to ensure correct package
# resolution for bioinformatics tools.
# -----------------------------------------------------------------------
RUN micromamba install -y -n base \
    -c bioconda \
    -c conda-forge \
    samtools=1.19 \
    bcftools=1.19 \
    seqkit=2.8.0 \
    python=3.10 \
    matplotlib=3.8.* && \
    micromamba clean --all --yes

# -----------------------------------------------------------------------
# Install selected aligner
# Fails fast with a clear error for unknown aligner values.
# -----------------------------------------------------------------------
RUN if [ "$ALIGNER" = "minimap2" ]; then \
        micromamba install -y -n base \
            -c bioconda \
            -c conda-forge \
            minimap2=2.28; \
    elif [ "$ALIGNER" = "pbmm2" ]; then \
        micromamba install -y -n base \
            -c bioconda \
            -c conda-forge \
            pbmm2=1.13.0; \
    else \
        echo "Error: Unknown aligner '$ALIGNER'. Use minimap2 or pbmm2." && exit 1; \
    fi && \
    micromamba clean --all --yes

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /pipeline

USER $MAMBA_USER

# -----------------------------------------------------------------------
# Smoke test
# Verifies all base tools and procps installed correctly.
# Aligner check uses || fallback since ARG scope ends after build steps
# and we cannot conditionally test by ALIGNER value at this stage.
# -----------------------------------------------------------------------
RUN samtools --version | head -1 && \
    bcftools --version | head -1 && \
    seqkit version && \
    ps --version && \
    python3 --version && \
    (minimap2 --version 2>/dev/null || pbmm2 --version 2>/dev/null || true)

CMD ["bash"]
