#!/bin/bash

if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    j=${SLURM_ARRAY_TASK_ID}
    d=$1
    SL_version=$2
else
    j=$1
    d=$2
    SL_version=$3
fi

Rscript END_senspar_cluster.R ${j} ${d} ${SL_version}
