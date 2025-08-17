# ![h3abionet/chipimputation](docs/images/h3abionet_logo.png)

[![GitHub Actions CI Status](https://github.com/h3abionet/chipimputation/workflows/nf-core%20CI/badge.svg)](https://github.com/h3abionet/chipimputation/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/h3abionet/chipimputation/workflows/nf-core%20linting/badge.svg)](https://github.com/h3abionet/chipimputation/actions?query=workflow%3A%22nf-core+linting%22)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/h3abionet/chipimputation)

## Introduction

**h3abionet/chipimputation** is a bioinformatics pipeline for genotype imputation and quality control.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation follows the [nf-core](https://nf-co.re) guidelines for best practices.

## Pipeline summary

1. Quality Control ([`BCFtools`](http://samtools.github.io/bcftools/))
   - Check and validate input VCF files
   - Remove duplicate variants
   - Split multi-allelic variants
   - Filter by minor allele count
   - Site missingness filtering

2. Phasing ([`Eagle2`](https://alkesgroup.broadinstitute.org/Eagle/))
   - Reference-based phasing using Eagle2
   - Chunking for parallel processing

3. Imputation ([`Minimac4`](https://genome.sph.umich.edu/wiki/Minimac4))
   - Genotype imputation using reference panels
   - Support for multiple reference panels

4. Post-imputation QC and Reporting
   - Filter by imputation quality score
   - Generate accuracy metrics
   - Create comprehensive QC plots
   - MultiQC report generation

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run h3abionet/chipimputation -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run h3abionet/chipimputation \
       -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
       --input samplesheet.csv \
       --outdir <OUTDIR>
   ```

## Documentation

The h3abionet/chipimputation pipeline comes with documentation about the pipeline [usage](https://github.com/h3abionet/chipimputation/blob/master/docs/usage.md), [parameters](https://github.com/h3abionet/chipimputation/blob/master/docs/parameters.md) and [output](https://github.com/h3abionet/chipimputation/blob/master/docs/output.md).

## Credits

h3abionet/chipimputation was originally written by:

- Mamana Mbiyavanga
- Gerrit Botha
- Eugene de Beste

We thank the following people for their extensive assistance in the development of this pipeline:

- The H3ABioNet Consortium
- The H3Africa Initiative
- nf-core community

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#chipimputation` channel](https://h3africa.slack.com/channels/chipimputation) (you can join with [this invite](https://h3africa.slack.com/signup)).

## Citations

If you use h3abionet/chipimputation for your analysis, please cite it using the following doi: [10.1038/s41467-018-05188-3](https://doi.org/10.1038/s41467-018-05188-3)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).