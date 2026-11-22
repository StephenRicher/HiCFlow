FROM condaforge/miniforge3:24.7.1-2

ARG NEXTFLOW_VERSION=24.10.4

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    NXF_HOME=/opt/nextflow \
    NXF_VER=${NEXTFLOW_VERSION} \
    PATH=/opt/conda/envs/histone-peaks/bin:/opt/conda/bin:${PATH}

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        gawk \
        gzip \
        openjdk-17-jre-headless \
        procps \
        tini \
    && rm -rf /var/lib/apt/lists/*

RUN mamba create \
        --yes \
        --name histone-peaks \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        python=3.11 \
        bedtools \
        bowtie2 \
        fastp \
        macs2=2.2.9.1 \
        pysradb \
        samtools \
        sra-tools \
        deepTools \
    && conda clean --all --yes

RUN cd /usr/local/bin \
    && curl \
        --fail \
        --location \
        --retry 10 \
        --retry-delay 10 \
        --silent \
        --show-error \
        https://get.nextflow.io \
        --output nextflow \
    && chmod +x nextflow \
    && nextflow -version

WORKDIR /workflow

COPY main.nf samples.tsv /workflow/

RUN python --version \
    && bedtools --version \
    && bowtie2 --version \
    && fastp --version \
    && macs2 --version \
    && pysradb --version \
    && samtools --version | head -n 1 \
    && fasterq-dump --version

ENTRYPOINT ["/usr/bin/tini", "--", "nextflow"]

CMD ["run", "/workflow/main.nf"]
