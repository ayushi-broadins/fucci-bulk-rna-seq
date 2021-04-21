#!/bin/bash

source /broad/software/scripts/useuse
reuse -q Anaconda
source activate /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/tools/fucciRNAenv
cd /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/scripts/

#build star index
#unzip the files
#gunzip "${1}.gz"
#gunzip "${2}.gz"


#create the star index
STAR --runMode genomeGenerate \
--genomeDir $3 \
--genomeFastaFiles $1 \
--sjdbGTFfile $2 \
--runThreadN 1
