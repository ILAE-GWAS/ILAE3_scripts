#!/bin/bash

conda activate ldsc

export APP_DIR=" "
export REF_DIR=" "
export EPI_DIR=" "
export SUMSTATS_DIR=" "
export OUT_DIR=" "
export TRAIT1=" "
export TRAIT2=" "

${APP_DIR}/ldsc.py \
--rg ${SUMSTATS_DIR}/${TRAIT1}.sumstats.gz,\
${SUMSTATS_DIR}/${TRAIT2}.sumstats.gz \
--ref-ld-chr ${REF_DIR}/ \
--w-ld-chr ${REF_DIR}/ \
--out ${OUT_DIR}/results \
