# HiCflow

## Comprehensive bioinformatics analysis pipeline for processing raw HiC read data to publication read HiC maps.

HiCflow aims to provide an accessible and user-friendly experience to analyse HiC data using a wide range of published tools.
The pipeline utilises the workflow management system Snakemake and automatically handles installation of all required software with no user input. HiCflow can also be easily scaled to work in cluster environments. Current software utilised by HiCflow include:

 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - A quality control tool for high throughput sequence data.
 * [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) - A tool to screen for species composition in FASTQ sequences.
 * [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - A tool to remove adapter sequences, primers, poly-A tails and others from high-throughput sequencing reads.
 * [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) - A tool for mapping and performing quality control on Hi-C data.
 * [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/) - A set of tools for building, normalising and processing HiC matrices.
 * [OnTAD](https://github.com/anlin00007/OnTAD) - An optimised nested TAD caller for identifying hierarchical TADs in HiC data.
 * [HiCRep](https://genome.cshlp.org/content/early/2017/08/30/gr.220640.117) - A tool for assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient.
 * [HiCcompare](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2288-x) - A tool joint normalisation and comparison of HI-C datasets
 * [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) - A tool for plotting customisable, publication ready genome tracks including HiC maps.
 * [MultiQC](https://multiqc.info/) - Aggregate results from bioinformatics analyses across many samples into a single report.

## Table of contents

  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Usage](#usage)
  * [Example output](#example-output)
     * [HiC track](#hic-track)
     * [HiCcompare track](#hiccompare-track)
  * [Quality Control](#quality-control)
     * [HiCRep](#hicrep)
     * [MultiQC report](#multiqc-report)
     * [Other QC metrics](#custom-qc-metrics)
  * [References](#references)

## Installation

HiCflow works with python >=3.6 and requires [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The HiCflow repository can be downloaded from github as follows:

```bash
$ git clone https://github.com/StephenRicher/HiCFlow.git
```

## Configuring HiCFlow

The HiCFlow pipeline is fully controlled through a single configuration file that describes parameter settings and paths to relevant files in the system.
HiCFlow is bundled with a fully configured small Hi-C dataset ([Wang et al., 2018](https://www.nature.com/articles/s41467-017-02526-9)) to test and serve as a template for configuring other datasets.
The configuration file for this example dataset is shown below and can be found at `example/config/config.yaml`.

**Note:** If relative file paths are provided in the configuration file then these are **relative to the working directory**.
The working directory itself (defined by workdir) is relative to the directory snakemake is executed in.
If not set, the working directory defaults to the directory containg the Snakefile.
Relative paths can be confusing, they are used here to ensure the example dataset works for all users.
If in doubt simply provide absolute paths.

```bash
# Specify output directory - either absolute path or relative to Snakefile.
# If using relative paths for subsequent files, these should be relative to
# this workding directory.
workdir: example/analysis/

# CSV file with cell type, experimental group, replicate number,
# read (forward/reverse) and path each FASTQ files.
data:  ../config/samples.csv

# Bed file of genomic regions to perform HiC analysis.
# These may be whole chromosomes for normal HiC or specific capture regions
# for region capture HiC.
regions: ../config/regions.bed

# FASTA references to align data. Must specify a reference for each cell type
# defined in config['data'].
genome :
    S2Rplus : ../genome/BDGP6.28.fa.gz

# Set True to perform phasing and haplotype assembly pipeline.
phase : True

# Phased VCF file for allele specific analysis. Must specify a VCF for each
# cell type defined in config['data']. If not set then run normal HiC mode.
# The HiCFlow phasing pipeline (see above) outputs a phased VCF for each cell
# type which is valid input here.
phased_vcf:
    #S2Rplus : ../analysis/phasedVCFs/S2Rplus-phased.vcf

# List of binsizes to analyse HiC data at different resolutions.
# The first binsize defines the base resolution, all subsequence bin sizes
# must be whole divisible by the base bin size e.g. [1000, 1500] is invalid.
binsize : [1000, 3000]

# Parameters for cutadapt - see https://cutadapt.readthedocs.io/en/stable/guide.html
cutadapt:
    forwardAdapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    reverseAdapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    overlap: 6
    errorRate: 0.1
    minimumLength: 20
    qualityCutoff: 20
    GCcontent: 43


# List of restriction sequence in order of protocol usage. Cut site is denoted
# using the '^' symbol. Ensure restriction enzyme names are given as strings.
restrictionSeqs:
    'DpnII' : '^GATC'

HiCParams:
    minDistance: 300
    maxLibraryInsertSize: 1000
    removeSelfLigation: True
    keepSelfCircles: False
    skipDuplicationCheck: False
    nofill: False

# Bigwig tracks for plotting below HiC plots.
bigWig :
    CP190 : ../genome/CP190-dm6.bw    # GSM762836
    Beaf-32 : ../genome/Beaf32-dm6.bw # GSM762845
    CTCF : ../genome/CTCF-dm6.bw      # GSM1535983
# BED tracks for plotting with HiC plots.
bed :
    Genes : ../genome/BDGP6.28.99.genes.bed

# BED file for creating plots of additional viewpoints in addition to those
# defined config['protocol']['regions'].
plot_coordinates: ../config/plot_coordinates.bed

# Matplotlib colour map for HiC plots
colourmap : 'Purples'

HiCcompare:
    fdr: 0.1 # FDR threshold for significant differential interactions.
    logFC: 0 # Fold-change threshold for significant differential interactions.
    multi: True # Run multiHiCcompare when replicates are available.

compareMatrices:
    vMin: -1.96 # Mimimum logFC value for colour scale.
    vMax: 1.96  # Maximum logFC value for colour scale.
    size: 3     # Size of median filter to denoise comparison matrix.

# Output a BAM file containing only valid HiC read pairs within defined regions.
createValidBam: False

# Configuration file for customising multiQC output report.
multiQCconfig : ../config/multiqc_config.yaml

# Configuration file of paths to genome indexes for FastQ Screen.
# See template in example/config/fastq_screen.config
fastq_screen :

```

## Usage

Once Snakemake is installed the example dataset can be processed using the following command.
This command should be run from the HiCFlow base directory, containing the Snakefile.

```bash
$ snakemake --use-conda --cores 4 --configfile example/config/config.yaml
```

This command will first install all of the relevant Conda environments within the defined working directory (`example/analysis/`).
This may take some time.
The pipeline should then run to completion producing the same figures as shown in the example output below.
Alteratively, you may also want to install the Conda environments in a custom directory.
This is useful when you are running multiple independent analyses and do not wish to repeatedly install the same Conda environments.

```bash
$ snakemake --use-conda --conda-prefix /path/envs/ --cores 4 --configfile example/config/config.yaml
```

### Cluster Execution
All Snakemake based pipelines, including HiCFlow, are compatible with cluster environments.
Consult the official Snakemake documentation [here](https://snakemake.readthedocs.io/en/v5.25.0/executing/cli.html#profiles) to learn more about running HiCFlow on your particular cluster environment.


## Example output

### HiC track

HiCflow utilises pyGenomeTracks to plot annotated HiC tracks with nested TAD domains, loops and TAD insulation scores. ChIP data and orientation of CTCF sites can also be provided.
![HiC plot example](./README_files/AS-chr3L-3L_5500000_6000000-3000.png)

### HiCcompare track

HiCflow uses HiCcompare to produce joint normalised Z-score subtraction matrices between pairs of samples. Statistically significant changes are highlighted as points on the subtraction matrix.
![HiCcompare example](./README_files/G1S-vs-AS-chr3L-3L_5500000_6000000-3000-logFC.png)

## Quality Control

### MultiQC report

HiCflow utilises MultiQC to aggregate the QC and metric report across all samples and all compatible tools used in the pipeline. An example MultiQC report produced by HiCflow is shown [here](./README_files/multiqc_report.html).  

### HiCRep

HiCflow uses HiCRep to assess sample-reproducibility by calculating the stratum-adjusted correlation coefficient between all pairwise samples.
![HiCRep example](./README_files/chr3L-1000-hicrep.png)

### Other QC Metrics

#### Insert Size Distribution
![InsertSize](./README_files/insert_size_frequency.png)

#### Ditag Length
![Ditag Length](./README_files/ditag_length.png)

###### References
Qi Wang, Qiu Sun, Daniel M. Czajkowsky, and Zhifeng Shao. Sub-kb Hi-C in D.
melanogaster reveals conserved characteristics of TADs between insect and mammalian
cells. Nature Communications, 2018. ISSN 20411723. doi: 10.1038/s41467-017-02526-9.
