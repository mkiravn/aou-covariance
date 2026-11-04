#!/usr/bin/env bash
# Local (macOS) setup for running the pipeline on simulated data.
# The Workbench has its own plink/plink2 under ~/bin; this is for a laptop.
set -euo pipefail
cd "$(dirname "$0")/../.."

python3 -m venv .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q numpy pandas scipy scikit-learn matplotlib jupyter papermill ipykernel nbstripout
.venv/bin/pip install -q -e .
# Python 3.13+ skips .pth files macOS has flagged hidden, which breaks the editable install
[ "$(uname -s)" = Darwin ] && chflags nohidden .venv/lib/python3.*/site-packages/*.pth
.venv/bin/nbstripout --install

Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  need <- c("bigsnpr", "bigutilsr", "data.table", "dplyr", "jsonlite", "e1071")
  need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org")'

mkdir -p .tools/bin
case "$(uname -m)" in
  arm64) P2=plink2_mac_arm64_20261024.zip ;;
  *)     P2=plink2_mac_20261024.zip ;;
esac
tmp=$(mktemp -d)
curl -sL -o "$tmp/p2.zip" "https://s3.amazonaws.com/plink2-assets/alpha7/$P2"
curl -sL -o "$tmp/p1.zip" "https://s3.amazonaws.com/plink1-assets/plink_mac_20260927.zip"
unzip -o -q "$tmp/p2.zip" plink2 -d .tools/bin
unzip -o -q "$tmp/p1.zip" plink -d .tools/bin
rm -rf "$tmp"
.tools/bin/plink2 --version
.tools/bin/plink --version

make -C src
echo "setup done"
