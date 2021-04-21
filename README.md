Bulk RNA-seq analysis pipeline with Snakemake for Broad Institute UGER cluster
======================================================================================

[Snakemake][snakemake] is a Pythonic workflow description language, that is 
easily configurable to run in all sorts of environments. Since version 4.1, 
Snakemake contains a feature called 'profiles', for easy exchange of 
configuration presets for running in a certain environment. This repository 
contains a snakemake bulk RNA-seq analysis pipeline to run on the Broad's UGER cluster.

<img src="https://github.com/ayushi-broadins/fucci-bulk-rna-seq/blob/13ce2946a9efbee264ba612822ddb99d465fbfa8/bulk_rna_seq_pipeline.png" >


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

1. To install RStudio and R, please follow the instructions [here][hr].
[hr]: https://uvastatlab.github.io/phdplus/installR.html
2. Open RStudio and install DESeq2 using the instructions below:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```


Using the pipeline
------------------------
We're ready to go! To start the snalysis:
1. Update the config file parameters/config.yaml to ensure it has the right paths 
   and sample names.
3. Connect to the Broad login host

```bash
ssh login

# tmux command will keep your code running even if you disconnect
# To disconnect from the session, press CTRL+b, release both keys 
# and then press d. On the original login shell, type tmux a to reconnect 
# to your session, tmux ls to list all sessions, and tmux a -t [number] to 
# connect to session [number].
tmux
cd scripts/
use Anaconda
source activate ../tools/snakemake
# The below command will run all the steps in the pipeline
# up until the alignment by STAR.
snakemake --profile broad-uger --cluster-config cluster.json

```
3. Open R Studio 
