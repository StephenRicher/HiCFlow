# Human Capture Hi-C analysis

Configuration for analysing the human chromosome 1 Capture Hi-C dataset using GRCh37/hg19.

## Required files

Large and analysis-specific files are not committed to Git. Add them under:

```text
analysis/
├── data/       # Paired FASTQ files
├── genome/     # hg19 FASTA and matching GFF3 annotation
└── tracks/     # BigWig tracks
```

Update `samples.tsv` and `analysis.config` with the real filenames.

All FASTA, GFF3, BED, and BigWig files must use GRCh37/hg19 coordinates and matching chromosome names. This configuration expects `chr1`, not `1`.

## Configuration

The reusable configuration is:

```text
analysis/
├── analysis.config
└── config/
    ├── samples.tsv
    ├── regions.bed
    ├── plot_coordinates.bed
    ├── viewpoints.bed
    ├── vlines.bed
    └── README.md
```

The exact Capture Hi-C interval is:

```text
chr1	66299999	70699997	chr1_capture
```

Configured plot windows cover the full capture interval and focused locus views. The viewpoints represent the `GNG12-AS1` TSS, B1 boundary, and `DIRAS3` deletion intervals.

## Download the reference files

This analysis uses GRCh37/hg19 with UCSC chromosome names such as `chr1`.

From the repository root:

```bash
mkdir -p analysis/genome

curl \
    --fail \
    --location \
    --retry 5 \
    --output analysis/genome/hg19.fa.gz \
    https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/latest/hg19.fa.gz

curl \
    --fail \
    --location \
    --retry 5 \
    --output analysis/genome/Homo_sapiens.GRCh37.75.gff3.gz \
    https://ftp.ensembl.org/pub/grch37/release-75/gff3/homo_sapiens/Homo_sapiens.GRCh37.75.gff3.gz
```

The UCSC FASTA uses chromosome names such as `chr1`, but the Ensembl GFF3 uses names such as `1`. Create a chromosome-prefixed annotation:

```bash
gzip -dc analysis/genome/Homo_sapiens.GRCh37.75.gff3.gz \
    | awk 'BEGIN {OFS="\t"}
        /^#/ {
            print
            next
        }
        $1 ~ /^[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "MT" {
            $1 = "chr" $1
        }
        {
            print
        }' \
    | gzip \
    > analysis/genome/Homo_sapiens.GRCh37.75.chr.gff3.gz
```

Verify that both files use `chr1`:

```bash
gzip -dc analysis/genome/hg19.fa.gz \
    | grep -m 1 '^>chr1'

gzip -dc analysis/genome/Homo_sapiens.GRCh37.75.chr.gff3.gz \
    | awk '$1 == "chr1" {print; exit}'
```

Every row in `samples.tsv` must use `Human` in the `cell_type` column so it matches the genome configuration key.

## Build

From the repository root:

```bash
docker build \
    --platform linux/amd64 \
    --tag hicflow-nextflow:latest \
    .
```

Configuration and input-file changes do not require rebuilding the image.

## Run

```bash
docker run --rm \
    --platform linux/amd64 \
    --volume "$PWD/analysis:/analysis" \
    --workdir /analysis \
    hicflow-nextflow:latest \
    run /opt/hicflow/main.nf \
    -profile container \
    -c /analysis/analysis.config \
    -work-dir /analysis/work \
    --outdir /analysis/results \
    -resume
```

Remove `-resume` for the first run if preferred. Keep `analysis/work/` to resume an interrupted run.

## Outputs

Results are written under:

```text
analysis/results/
├── bam/
├── loops/
├── matrix/
├── pipeline_info/
├── plots/
├── qc/
├── tads/
└── viewpoints/
```

Subtraction direction is explicit in filenames. For example:

```text
group_2-minus-group_1
```

means `group_2 - group_1`. Positive values are higher in `group_2`; negative values are higher in `group_1`.

## Clean rerun

Delete cached work and outputs only when a completely fresh run is required:

```bash
rm -rf analysis/work analysis/results
```

Then rerun the command without `-resume`.

Only `analysis/analysis.config` and `analysis/config/` should be committed. FASTQs, references, tracks, work files, logs, and results should remain ignored.
