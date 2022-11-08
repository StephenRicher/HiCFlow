# HiCflow

## Comprehensive bioinformatics analysis pipeline for processing raw HiC read data to publication read HiC maps.

HiCflow aims to provide an accessible and user-friendly experience to analyse HiC data using a wide range of published tools.
The pipeline utilises the workflow management system Snakemake and automatically handles the installation of all required software with no user input. HiCflow can also be easily scaled to work in cluster environments. Current software utilised by HiCflow includes:

 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - A quality control tool for high throughput sequence data.
 * [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) - A tool to screen for species composition in FASTQ sequences.
 * [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - A tool to remove adapter sequences, primers, poly-A tails and others from high-throughput sequencing reads.
 * [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) - A tool for mapping and performing quality control on HiC data.
 * [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/) - A set of tools for building, normalising and processing HiC matrices.
 * [OnTAD](https://github.com/anlin00007/OnTAD) - An optimised nested TAD caller for identifying hierarchical TADs in HiC data.
 * [HiCRep](https://genome.cshlp.org/content/early/2017/08/30/gr.220640.117) - A tool for assessing the reproducibility of HiC data using a stratum-adjusted correlation coefficient.
 * [HiCcompare](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2288-x) - A tool for joint normalisation and comparison of HI-C datasets
 * [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) - A tool for plotting customisable, publication-ready genome tracks including HiC maps.
 * [MultiQC](https://multiqc.info/) - Aggregate results from bioinformatics analyses across many samples into a single report.

## Table of contents

  * [Installation](#installation)
  * [Configuration](#configuration)
    * [Example Configurations](#example-configurations)
  * [Usage](#usage)
  * [Example output](#example-output)
     * [HiC track](#hic-track)
     * [HiCcompare track](#hiccompare-track)
     * [Viewpoints](#viewpoints)
  * [Quality Control](#quality-control)
     * [HiCRep](#hicrep)
     * [MultiQC report](#multiqc-report)
     * [Other QC metrics](#custom-qc-metrics)
  * [References](#references)

## Installation

HiCflow works with python >=3.6 and requires [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The HiCflow repository can be downloaded from GitHub as follows:

```bash
git clone https://github.com/StephenRicher/HiCFlow.git
```

## Configuring HiCFlow

The HiCFlow pipeline is fully controlled through a single configuration file that describes parameter settings and paths to relevant files in the system.
HiCFlow is bundled with a fully configured small HiC dataset ([Wang et al., 2018](https://www.nature.com/articles/s41467-017-02526-9)) to test and serve as a template for configuring other datasets.
The configuration file for the example dataset is provided at at [here](example/config/config.yaml).

**Note:** If relative file paths are provided in the configuration file, then these are **relative to the working directory**.
The working directory itself (defined by workdir) is relative to the directory ``snakemake`` is executed.
If not set, the working directory defaults to the directory containing the Snakefile.
Relative paths can be confusing; they are used here to ensure the example dataset works for all users.
If in doubt, simply provide absolute paths.

### Example Configurations

* [Typical HiC Analysis](example/config/config-HiC+CallVariant+Phase.yaml)
  * Run standard HiC workflow.
* [HiC Analysis + Variant Calling + Haplotype Assembly](example/config/config-HiC+CallVariant+Phase.yaml)
  * Run standard HiC workflow and full variant calling and haplotype assembly pipeline.
  * Phased VCF output compatible with ASHiC workflow.
* [HiC Analysis + Haplotype Assembly](example/config/config-HiC+Phase.yaml)
  * Run standard HiC workflow and haplotype assembly pipeline.
  * Requires a set of pre-called variants.
  * **If high quality calls from WGS data are available then we recommended using these rather than performing variant calling from the HiC data using HiCFlow.**
* [Allele Specific HiC](example/config/config-ASHiC.yaml)
  * Perform allele-specific HiC workflow.
  * Requires a set of phased variants, either from HiCFlow or another source.

## Usage

Once Snakemake is installed, the example dataset can be processed using the following command.
This command should be run from the HiCFlow base directory containing the Snakefile.

```bash
snakemake --use-conda --cores 4 --configfile example/config/config.yaml
```

This command will first install all relevant Conda environments within the defined working directory (`example/analysis/`); this may take some time.
The pipeline should then run to completion producing the exact figures as shown in the example output below.
Alternatively, you may also want to install the Conda environments in a custom directory.
A custom directory is helpful if you perform multiple independent analyses and do not want to install the same Conda environments repeatedly.

```bash
snakemake --use-conda --conda-prefix /path/envs/ --cores 4 --configfile example/config/config.yaml
```

### Cluster Execution
All Snakemake-based pipelines, including HiCFlow, are compatible with cluster environments.
Consult the official Snakemake documentation [here](https://snakemake.readthedocs.io/en/v5.25.0/executing/cli.html#profiles) to learn more about running HiCFlow on your particular cluster environment.


## Example output

### HiC track

HiCflow utilises pyGenomeTracks to plot annotated HiC tracks with nested TAD domains, loops and TAD insulation scores. In addition, custom BED and Bedgraph files can be provided through the configuration file.
![HiC plot example](./README_files/AS-chr3L-3L_5500000_6000000-3000-custom-full-fm.png)

### HiCcompare track

HiCflow uses HiCcompare to produce joint normalised log fold-change subtraction matrices between pairs of samples.
![HiCcompare example](./README_files/G1S-vs-AS-chr3L-3L_5500000_6000000-3000-logFC-full-fm.png)

### Viewpoints

HiCFlow can also plot custom viewpoints of specific regions. Viewpoint regions must be provided as a BED file in the configuration file under plotParams -> viewpoints. The below example compares two samples using between-sample normalised contact frequencies provided by HiCcompare.
![Viewpoint example](./README_files/G1S-vs-AS-chr3L-3L_5740000_5750000-3000-viewpoint-full.svg)

## Quality Control

### MultiQC report

HiCflow utilises MultiQC to aggregate the QC and metric report across all samples and all compatible tools used in the pipeline. An example MultiQC report produced by HiCflow is shown [here](./README_files/multiqc_report.html).  

### HiCRep

HiCflow uses HiCRep to assess sample reproducibility by calculating the stratum-adjusted correlation coefficient between all pairwise samples.
![HiCRep example](./README_files/chr3L-1000-hicrep-full.svg)

### Other QC Metrics

#### Insert Size Distribution
![InsertSize](./README_files/insertSizeFrequency.svg)

#### Ditag Length
![Ditag Length](./README_files/ditagLength.svg)

## References
Qi Wang, Qiu Sun, Daniel M. Czajkowsky, and Zhifeng Shao. Sub-kb Hi-C in D.
melanogaster reveals conserved characteristics of TADs between insect and mammalian
cells. Nature Communications, 2018. ISSN 20411723. doi: 10.1038/s41467-017-02526-9.
