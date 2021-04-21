#!/bin/bash

source /broad/software/scripts/useuse
reuse -q UGER
reuse -q R-3.5
cd /chembio/datasets/csdev/AA/betacellproliferation/rnaseq_fucci/scripts/

#run the r script
Rscript r/parse.read.counts.star.R $1 $2 $3