#!/usr/bin/env bash
# Submit one dsub task array to compute all N GRM shards, modeled on the
# aou_dsub + --tasks pattern from asgarilab/SDoH-GeneticRisk-Biobank-MCA.
# Run this from a Jupyter/RStudio cell on the AoU Researcher Workbench,
# where `aou_dsub` (~/aou_dsub.bash) and $WORKSPACE_BUCKET are already set up.
#
# Usage: bash submit_grm_shards.sh <bfile_gcs_dir> <n_shards> <out_prefix>
set -euo pipefail

BFILE_DIR="$1"      # e.g. gs://fc-aou-datasets-controlled/.../plink_v7.1
N_SHARDS="$2"
OUT_PREFIX="$3"     # e.g. grm_ceugbr

PARAM_FILE="grm_shards_param.tsv"
{
    echo -e "--env SHARD_K"
    for k in $(seq 1 "${N_SHARDS}"); do
        echo -e "${k}"
    done
} > "${PARAM_FILE}"

cat > script.sh <<'EOS'
#!/bin/bash
set -euo pipefail
plink --bfile "${bfile_dir}/arrays" \
    --make-grm-bin --parallel "${SHARD_K}" "${N_SHARDS}" \
    --nonfounders \
    --out "${outpath}/${OUT_PREFIX}"
EOS

gsutil cp script.sh "${WORKSPACE_BUCKET}/data/scripts/bash/grm_shards/submit_grm_shards.sh"
gsutil cp "${PARAM_FILE}" "${WORKSPACE_BUCKET}/data/meta/grm_shards/${PARAM_FILE}"

source ~/aou_dsub.bash

aou_dsub \
    --image biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    --min-ram 16 --min-cores 8 \
    --input-recursive bfile_dir="${BFILE_DIR}" \
    --output-recursive outpath="${WORKSPACE_BUCKET}/data/grm_shards" \
    --env N_SHARDS="${N_SHARDS}" --env OUT_PREFIX="${OUT_PREFIX}" \
    --logging "${WORKSPACE_BUCKET}/data/logs/{job-name}/{job-id}-{task-id}.log" \
    --tasks "${WORKSPACE_BUCKET}/data/meta/grm_shards/${PARAM_FILE}" \
    --script "${WORKSPACE_BUCKET}/data/scripts/bash/grm_shards/submit_grm_shards.sh"

echo "Submitted ${N_SHARDS} shard jobs. Monitor with dstat as in GRM-pairs/aou/batch/ examples."
