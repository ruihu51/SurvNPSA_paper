#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

usage() {
    echo "Usage:"
    echo "  bash submit_sim.sh <j> <n> [gam|superlearner]"
    echo "  sbatch --array=501-1500%100 submit_sim.sh <n> [gam|superlearner]"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    j="${SLURM_ARRAY_TASK_ID}"
    N="${1:-}"
    method="${2:-gam}"
else
    j="${1:-}"
    N="${2:-}"
    method="${3:-gam}"
fi

if [[ -z "$j" || -z "$N" ]]; then
    usage >&2
    exit 1
fi

if [[ "$method" != "gam" && "$method" != "superlearner" ]]; then
    echo "Error: method must be either 'gam' or 'superlearner'." >&2
    usage >&2
    exit 1
fi

Rscript sim_main.R "$j" "$N" "$method"
