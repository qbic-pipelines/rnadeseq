# qbic-pipelines/rnadeseq

**Downstream differential gene expression analysis with DESeq2 package**.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/qbicpipelines/rnadeseq.svg)](https://hub.docker.com/r/qbicpipelines/rnadeseq)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4024184.svg)](https://doi.org/10.5281/zenodo.4024184)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Documentation

The qbic-pipelines/rnadeseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

qbic-pipelines/rnadeseq was originally written by Gisela Gabernet (@ggabernet) and Silvia Morini (@silviamorins), at QBiC.

The pipeline structure is based on the template by the `nf-core` project. For more information, please check out the [nf-core website](https://nf-co.re/).
