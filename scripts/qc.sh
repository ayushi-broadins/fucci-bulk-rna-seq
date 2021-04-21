#!/bin/bash

source /broad/software/scripts/useuse
reuse -q Anaconda
source activate /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/tools/fucciRNAenv
cd /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/scripts/

#run fastqc for the fastq files
fastqc --outdir $1 $2