FROM condaforge/miniforge3:24.7.1-2

ARG NEXTFLOW_VERSION=24.10.4

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    NXF_HOME=/opt/nextflow \
    NXF_VER=${NEXTFLOW_VERSION} \
    PATH=/opt/conda/envs/corigami-preprocess/bin:/opt/conda/bin:${PATH}

RUN sed -i 's|http://archive.ubuntu.com|https://archive.ubuntu.com|g' /etc/apt/sources.list \
    && sed -i 's|http://security.ubuntu.com|https://security.ubuntu.com|g' /etc/apt/sources.list \
    && apt-get \
        -o Acquire::Retries=10 \
        -o Acquire::https::Timeout=120 \
        update \
    && apt-get \
        -o Acquire::Retries=10 \
        -o Acquire::https::Timeout=120 \
        install \
        --yes \
        --no-install-recommends \
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
        --name corigami-preprocess \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        python=3.11 \
        bowtie2 \
        deeptools \
        fastp \
        pysradb \
        samtools \
        sra-tools \
    && conda clean --all --yes

RUN mamba create \
        --yes \
        --name corigami \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        python=3.9 \
        pip \
        "numpy<2" \
        pandas \
        matplotlib \
        omegaconf \
        tqdm \
        pybigwig=0.3.18 \
    && /opt/conda/envs/corigami/bin/python -m pip install \
        --index-url https://download.pytorch.org/whl/cpu \
        torch==2.2.2 \
    && /opt/conda/envs/corigami/bin/python -m pip install \
        --no-deps \
        corigami==0.0.1 \
    && /opt/conda/envs/corigami/bin/python -c \
        "import corigami, pyBigWig, torch; print(f'torch={torch.__version__}'); print(f'cuda={torch.cuda.is_available()}')" \
    && conda clean --all --yes \
    && rm -rf /root/.cache/pip

RUN cd /usr/local/bin \
    && curl \
        --fail \
        --location \
        --retry 10 \
        --retry-delay 10 \
        --retry-all-errors \
        --silent \
        --show-error \
        https://get.nextflow.io \
        --output nextflow \
    && chmod +x nextflow \
    && nextflow -version

WORKDIR /workflow

COPY \
    main.nf \
    nextflow.config \
    samples_*.tsv \
    deletions_*.tsv \
    /workflow/

COPY modules /workflow/modules

COPY scripts /workflow/scripts

RUN python --version \
    && bowtie2 --version \
    && fastp --version \
    && samtools --version | head -n 1 \
    && bamCoverage --version \
    && bigwigCompare --version \
    && pysradb --version \
    && fasterq-dump --version

RUN /opt/conda/envs/corigami/bin/corigami-predict --help >/dev/null \
    && /opt/conda/envs/corigami/bin/corigami-edit --help >/dev/null

ENTRYPOINT ["/usr/bin/tini", "--", "nextflow"]

CMD ["run", "/workflow/main.nf", "-c", "/workflow/nextflow.config"]
