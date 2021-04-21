Bulk RNA-seq analysis pipeline with Snakemake for Broad Institute UGER cluster
======================================================================================

[Snakemake][snakemake] is a Pythonic workflow description language, that is 
easily configurable to run in all sorts of environments. Since version 4.1, 
Snakemake contains a feature called 'profiles', for easy exchange of 
configuration presets for running in a certain environment. This repository 
contains a snakemake bulk RNA-seq analysis pipeline to run on the Broad's UGER cluster.

[snakemake]: https://snakemake.readthedocs.io/

Installation
------------

### Setting up snakemake profile for Broad Institute UGER cluster
Please follow the instructions on the Broad Institute [GitHub page][gp] to set up the 
snakemake profile.

[gp]: https://github.com/broadinstitute/snakemake-broad-uger/blob/master/README.md

### Preparing a conda environment
The recommended way to run the analysis in this repository is to setup a conda environment, 
where the package versions of the tools used can be controlled. For a windows machine, 
please follow the below installation instructions on the Broad cluster.

```bash
use Anaconda3

# Create new conda environment with the environment.yml file provided in this repository
dos2unix environment.yml
conda env create -f environment.yml

```
### Installing R Studio and DESeq2 package
[DESeq2][deseq2] is a R Bioconductor package that is used for differential expression analysis. 
This tool allows you to have more than two experimental groups and account for a second 
experimental factor. This tool takes as input a table of raw counts. 

[deseq2]: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

1. To install RStudio and R, please follow the instructions .[here].[hr].
.[hr]: https://uvastatlab.github.io/phdplus/installR.html
2. To install DESeq2, please follow the instructions below:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```


Using the pipeline
------------------------
We're ready to go! To start the snalysis:
