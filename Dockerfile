FROM condaforge/miniforge3:24.7.1-2

ARG NEXTFLOW_VERSION=24.10.4

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    NXF_HOME=/opt/nextflow \
    NXF_VER=${NEXTFLOW_VERSION} \
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

COPY environment.yml /tmp/environment.yml

RUN /opt/conda/bin/mamba env create \
        --file /tmp/environment.yml \
        --prefix /opt/conda/envs/hicflow \
    && /opt/conda/bin/conda clean --all --yes \
    && rm /tmp/environment.yml

WORKDIR /opt/hicflow

COPY . /opt/hicflow

RUN chmod +x \
        /opt/hicflow/workflow/scripts/hicup_truncater \
        /opt/hicflow/workflow/scripts/OnTAD \
        /opt/hicflow/workflow/scripts/CscoreTool1.1

RUN fastqc --version \
    && cutadapt --version \
    && bowtie2 --version \
    && samtools --version | head -n 1 \
    && hicBuildMatrix --version \
    && hicCorrectMatrix --version \
    && pyGenomeTracks --version \
    && multiqc --version \
    && Rscript -e "library(HiCcompare); library(hicrep)"

ENTRYPOINT ["/usr/bin/tini", "--", "nextflow"]

CMD ["run", "/opt/hicflow/main.nf", "-profile", "container"]
