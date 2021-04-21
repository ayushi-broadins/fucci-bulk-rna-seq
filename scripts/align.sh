#!/bin/bash

source /broad/software/scripts/useuse
reuse -q Anaconda
source activate /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/tools/fucciRNAenv
cd /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/scripts/

#align reads to genome using STAR aligner
STAR --runThreadN $1 \
--twopassMode Basic \
--quantMode GeneCounts \
--alignSJDBoverhangMin 1 \
--seedSearchStartLmax 30 \
--outFilterType BySJout \
--genomeDir $2 \
--readFilesIn $3 $4 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $5