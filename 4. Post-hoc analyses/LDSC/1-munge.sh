#!/bin/bash

# munge sumstats for LDSC genetic correlations

conda activate ldsc

export APP_DIR=" "
export TUT_DIR=" "
export SS_DIR=" "
export OUT_DIR=" "
export TRAIT=" "

## PREPARE SUMSTAT FILES ##

${APP_DIR}/munge_sumstats.py \
--sumstats ${SS_DIR}/${TRAIT}_sumstats.txt.gz \
--snp rsid \
--a1 EA \
--a2 NEA \
--p P-value \
--N-col Effective_N \
--out ${OUT_DIR}/${TRAIT} \
--merge-alleles ${TUT_DIR}/w_hm3.snplist
