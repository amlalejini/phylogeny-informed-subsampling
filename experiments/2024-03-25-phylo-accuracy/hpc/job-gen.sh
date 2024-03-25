#!/usr/bin/env bash

REPLICATES=10
EXP_SLUG=2024-03-25-phylo-accuracy
SEED_OFFSET=70000
JOB_TIME=48:00:00
JOB_MEM=16G
ACCOUNT=NONE
PROJECT_NAME=phylogeny-informed-subsampling

# SCRATCH_EXP_DIR=/mnt/scratch/lalejini/data-public/${PROJECT_NAME}
# REPO_DIR=/mnt/home/lalejini/devo_ws/${PROJECT_NAME}
# HOME_EXP_DIR=${REPO_DIR}/experiments
# TODO - Where to dump? Where to hold repo?

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR_SRC=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config
CONFIG_DIR=${JOB_DIR}/config

mkdir -p ${CONFIG_DIR}
cp -R ${CONFIG_DIR_SRC}/* ${CONFIG_DIR}/

python3 gen-sub.py --runs_per_subdir 500 --time_request ${JOB_TIME} --mem ${JOB_MEM} --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --repo_dir ${REPO_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}