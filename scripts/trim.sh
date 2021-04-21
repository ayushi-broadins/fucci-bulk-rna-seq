#!/bin/bash

source /broad/software/scripts/useuse
reuse -q Anaconda
source activate /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/tools/fucciRNAenv
cd /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/scripts/

#trim the illumina adapters from the fastq files
trimmomatic PE -phred33 $1 $2 $3 $4 $5 $6 ILLUMINACLIP:../tools/fucciRNAenv/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36