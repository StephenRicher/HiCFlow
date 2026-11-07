FROM condaforge/miniforge3:24.7.1-2

ARG NEXTFLOW_VERSION=24.10.4

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    NXF_HOME=/opt/nextflow \
    PATH=/opt/conda/envs/hicflow/bin:/opt/conda/bin:${PATH}

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        fontconfig \
        git \
        openjdk-17-jre-headless \
        procps \
        tini \
        unzip \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml

RUN /opt/conda/bin/mamba env create \
        --file /tmp/environment.yml \
        --prefix /opt/conda/envs/hicflow \
    && /opt/conda/bin/conda clean --all --yes \
    && rm /tmp/environment.yml

RUN curl --fail --location --silent --show-error \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow" \
        --output /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

WORKDIR /opt/hicflow

COPY . /opt/hicflow

RUN chmod +x /opt/hicflow/workflow/scripts/hicup_truncater

RUN fastqc --version \
    && cutadapt --version \
    && bowtie2 --version \
    && samtools --version | head -n 1 \
    && hicBuildMatrix --version \
    && hicCorrectMatrix --version \
    && pyGenomeTracks --version \
    && multiqc --version \
    && Rscript --version

ENTRYPOINT ["/usr/bin/tini", "--", "nextflow"]

CMD ["run", "/opt/hicflow/main.nf", "-profile", "container"]
