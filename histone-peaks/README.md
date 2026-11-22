# Histone and CTCF peak workflow

This standalone Nextflow workflow generates reproducible histone-mark and CTCF tracks for the human HiCFlow analysis.

It processes public MCF7 and MCF10A datasets for:

- H3K27ac
- H3K4me1
- CTCF

MCF10A is used as a non-transformed mammary epithelial comparator where equivalent public 1-7HB2 data are unavailable.

## Workflow

For each biological replicate, the workflow:

1. Resolves the GEO accession to its SRA run or runs.
2. Downloads the FASTQ files.
3. Downloads the UCSC hg19 reference.
4. Trims reads with fastp.
5. Aligns reads to hg19 with Bowtie2.
6. Retains primary alignments with MAPQ ≥30.
7. Removes duplicates with samtools.
8. Calls peaks with MACS2:
   - Broad peaks for H3K27ac and H3K4me1
   - Narrow peaks for CTCF
9. Retains peaks detected in both biological replicates.
10. Restricts the final peak and signal tracks to:

```text
chr1:66,299,999-70,699,997
```

The workflow also creates CPM-normalised BigWig tracks and averages the two replicates for each cell line and mark.

## Input datasets

The samples are defined in `samples.tsv`.

The workflow currently processes two replicates for each combination of:

```text
MCF7    H3K27ac
MCF7    H3K4me1
MCF7    CTCF
MCF10A  H3K27ac
MCF10A  H3K4me1
MCF10A  CTCF
```

## Build

Run from `HiCFlow/histone-peaks`:

```bash
docker build \
    --platform linux/amd64 \
    --tag hicflow-histone-peaks:latest \
    .
```

## Run

```bash
docker run --rm \
    --platform linux/amd64 \
    --volume "$PWD/results:/data" \
    --workdir /data \
    hicflow-histone-peaks:latest \
    run /workflow/main.nf \
    -work-dir /data/work \
    --outdir /data/results \
    -resume
```

The corresponding host directories are:

```text
results/work/      Nextflow work cache
results/results/   Published outputs
```

Keep the mount and work-directory arguments unchanged when using `-resume`, otherwise Nextflow will not find the existing cache.

## Final peak tracks

The main capture-region BED files are written to:

```text
results/results/peaks/consensus/
├── MCF7_CTCF_chr1_capture.hg19.bed
├── MCF7_H3K27ac_chr1_capture.hg19.bed
├── MCF7_H3K4me1_chr1_capture.hg19.bed
├── MCF10A_CTCF_chr1_capture.hg19.bed
├── MCF10A_H3K27ac_chr1_capture.hg19.bed
└── MCF10A_H3K4me1_chr1_capture.hg19.bed
```

Here, `consensus` means that a peak called independently in replicate 1 overlaps a peak called in replicate 2.

The full-genome consensus files are also retained:

```text
results/results/peaks/consensus/
├── MCF7_CTCF_consensus.hg19.bed
├── MCF7_H3K27ac_consensus.hg19.bed
├── MCF7_H3K4me1_consensus.hg19.bed
├── MCF10A_CTCF_consensus.hg19.bed
├── MCF10A_H3K27ac_consensus.hg19.bed
└── MCF10A_H3K4me1_consensus.hg19.bed
```

The `*_chr1_capture.hg19.bed` files are the tracks intended for the HiCFlow plots.

## Final signal tracks

Mean CPM-normalised BigWigs are written to:

```text
results/results/signal/consensus/
├── MCF7_CTCF_mean.CPM.bw
├── MCF7_H3K27ac_mean.CPM.bw
├── MCF7_H3K4me1_mean.CPM.bw
├── MCF10A_CTCF_mean.CPM.bw
├── MCF10A_H3K27ac_mean.CPM.bw
└── MCF10A_H3K4me1_mean.CPM.bw
```

These BigWigs contain signal only for the configured capture interval.

## Other outputs

```text
results/results/bam/
```

Filtered, coordinate-sorted BAMs and their indexes.

```text
results/results/metadata/
```

GEO-to-SRA mappings used to retrieve the raw reads.

```text
results/results/peaks/replicates/
```

Per-replicate MACS2 peak calls normalised to BED3.

```text
results/results/qc/fastp/
```

Per-replicate fastp HTML and JSON reports.

```text
results/results/qc/bowtie2/
```

Per-replicate Bowtie2 alignment logs.

```text
results/results/logs/macs2/
```

Per-replicate MACS2 logs.

```text
results/results/signal/replicates/
```

Per-replicate CPM-normalised BigWigs restricted to the capture interval.

## Copy tracks into the HiCFlow analysis

From the `HiCFlow` repository root:

```bash
mkdir -p analysis/tracks

cp histone-peaks/results/results/peaks/consensus/*_chr1_capture.hg19.bed \
    analysis/tracks/

cp histone-peaks/results/results/signal/consensus/*.CPM.bw \
    analysis/tracks/
```

The resulting files under `analysis/tracks/` should be:

```text
analysis/tracks/
├── MCF7_CTCF_chr1_capture.hg19.bed
├── MCF7_H3K27ac_chr1_capture.hg19.bed
├── MCF7_H3K4me1_chr1_capture.hg19.bed
├── MCF10A_CTCF_chr1_capture.hg19.bed
├── MCF10A_H3K27ac_chr1_capture.hg19.bed
├── MCF10A_H3K4me1_chr1_capture.hg19.bed
├── MCF7_CTCF_mean.CPM.bw
├── MCF7_H3K27ac_mean.CPM.bw
├── MCF7_H3K4me1_mean.CPM.bw
├── MCF10A_CTCF_mean.CPM.bw
├── MCF10A_H3K27ac_mean.CPM.bw
└── MCF10A_H3K4me1_mean.CPM.bw
```

## HiCFlow configuration

Configure the peak tracks as BED tracks:

```groovy
bedTracks = [
    'MCF7 CTCF peaks'       : '/analysis/tracks/MCF7_CTCF_chr1_capture.hg19.bed',
    'MCF7 H3K27ac peaks'    : '/analysis/tracks/MCF7_H3K27ac_chr1_capture.hg19.bed',
    'MCF7 H3K4me1 peaks'    : '/analysis/tracks/MCF7_H3K4me1_chr1_capture.hg19.bed',
    'MCF10A CTCF peaks'     : '/analysis/tracks/MCF10A_CTCF_chr1_capture.hg19.bed',
    'MCF10A H3K27ac peaks'  : '/analysis/tracks/MCF10A_H3K27ac_chr1_capture.hg19.bed',
    'MCF10A H3K4me1 peaks'  : '/analysis/tracks/MCF10A_H3K4me1_chr1_capture.hg19.bed'
]
```

Configure the mean signal tracks separately:

```groovy
bigWig = [
    'MCF7 CTCF signal'      : '/analysis/tracks/MCF7_CTCF_mean.CPM.bw',
    'MCF7 H3K27ac signal'   : '/analysis/tracks/MCF7_H3K27ac_mean.CPM.bw',
    'MCF7 H3K4me1 signal'   : '/analysis/tracks/MCF7_H3K4me1_mean.CPM.bw',
    'MCF10A CTCF signal'    : '/analysis/tracks/MCF10A_CTCF_mean.CPM.bw',
    'MCF10A H3K27ac signal' : '/analysis/tracks/MCF10A_H3K27ac_mean.CPM.bw',
    'MCF10A H3K4me1 signal' : '/analysis/tracks/MCF10A_H3K4me1_mean.CPM.bw'
]
```

The BED files provide reproducible peak calls. The BigWigs provide supporting enrichment profiles.

## Validation

Count the final BED tracks:

```bash
find results/results/peaks/consensus \
    -maxdepth 1 \
    -type f \
    -name '*_chr1_capture.hg19.bed' \
    -print \
    | sort
```

Expected: six files.

Count the final merged BigWigs:

```bash
find results/results/signal/consensus \
    -maxdepth 1 \
    -type f \
    -name '*_mean.CPM.bw' \
    -print \
    | sort
```

Expected: six files.

Check that every capture-region BED is non-empty:

```bash
for file in results/results/peaks/consensus/*_chr1_capture.hg19.bed; do
    test -s "$file" || {
        echo "Missing or empty peak track: $file" >&2
        exit 1
    }

    printf '%s\t%s peaks\n' \
        "$file" \
        "$(wc -l < "$file" | tr -d ' ')"
done
```

## Interpretation

The consensus BED files represent peaks reproducible across the two selected biological replicates.

MCF7 and MCF10A were generated in different experimental contexts, including CUT&Tag and ChIP-seq. The tracks are suitable for asking whether CTCF or enhancer-associated histone marks are present around candidate regulatory elements.

Do not treat the following as directly quantitative between cell lines:

- Peak counts
- Peak widths
- MACS2 scores
- Absolute BigWig heights

The defensible interpretation is qualitative evidence of reproducible CTCF occupancy or enhancer-associated chromatin within each cell-line context.

## Track generation and provenance
Public sequencing data were downloaded from the NCBI Gene Expression Omnibus on 18 November 2026.
Wild-type MCF7 H3K27ac and H3K4me1 CUT&Tag replicates were obtained from GSE201262; untreated 0-hour MCF10A-ER-Src H3K27ac and H3K4me1 ChIP-seq replicates were obtained from GSE100257.
Normoxic MCF7 CTCF ChIP-seq replicates GSM8152452 and GSM8152453 and MCF10A CTCF ChIP-seq replicates GSM2599084 and GSM2599085 were analysed alongside the histone-mark datasets.
GEO accessions were resolved to SRA runs using pysradb, and reads were retrieved with prefetch and fasterq-dump.
Reads were adapter- and quality-trimmed with fastp, aligned independently to the UCSC hg19/GRCh37 reference using Bowtie2 in sensitive mode, and filtered with SAMtools to retain primary alignments with mapping quality ≥30.
Alignments were coordinate-sorted and PCR duplicates removed with samtools markdup.
Peaks were called independently for each replicate using MACS2 2.2.9.1 without matched input controls.
CTCF was called as narrow peaks using --qvalue 0.01 --nomodel --shift 0 --extsize 200; H3K27ac and H3K4me1 were called as broad peaks using --broad --broad-cutoff 0.1 --qvalue 0.05 --nomodel --shift 0 --extsize 200.
Consensus tracks comprise replicate-1 peak intervals overlapping at least one independently called replicate-2 peak; retained intervals were reduced to BED3, sorted and merged with BEDTools.
Consensus peaks were restricted to the hg19 capture interval chr1:66,299,999-70,699,997.
CPM-normalised 25 bp BigWig signal tracks were generated from each filtered BAM using deepTools bamCoverage, restricted to the same interval, and the two replicates for each cell line and mark were averaged using bigwigCompare --operation mean.
MCF10A was used as a non-transformed mammary epithelial comparator where equivalent public 1-7HB2 chromatin-profile data were unavailable.
Because the MCF7 and MCF10A datasets were generated using different assays and experimental contexts, track comparisons were interpreted qualitatively in terms of reproducible peak presence and genomic position rather than absolute signal magnitude.
