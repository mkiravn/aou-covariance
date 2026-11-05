#!/usr/bin/env bash
# Workbench (Linux) setup: Python and R packages into the Jupyter environment, plink binaries
# into ~/bin if missing, grm_class_tool. Run from a terminal on the VM.
set -euo pipefail
cd "$(dirname "$0")/../.."

python3 -m pip install -q numpy pandas scipy scikit-learn matplotlib papermill ipykernel nbformat
python3 -m pip install -q -e .

Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  need <- c("bigsnpr", "bigutilsr", "data.table", "dplyr", "jsonlite", "e1071")
  need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org")'

mkdir -p ~/bin
tmp=$(mktemp -d)
if [ ! -x ~/bin/plink2 ]; then
  curl -fsSL -o "$tmp/p2.zip" https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20261024.zip
  unzip -o -q "$tmp/p2.zip" plink2 -d ~/bin
fi
if [ ! -x ~/bin/plink ]; then
  curl -fsSL -o "$tmp/p1.zip" https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip
  unzip -o -q "$tmp/p1.zip" plink -d ~/bin
fi
rm -rf "$tmp"
~/bin/plink2 --version
~/bin/plink --version

make -C src
command -v zstdcat >/dev/null || echo "warning: zstd missing (hm3_overlap needs it): sudo apt-get install -y zstd"
command -v gcloud >/dev/null || echo "warning: gcloud missing (aoucov sync needs it)"
df -h ~ | tail -1
echo "setup done"
