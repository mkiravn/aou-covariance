#!/usr/bin/env bash
# Workbench (Linux) setup: Python and R packages into the Jupyter environment, plink binaries
# into ~/bin if missing, grm_class_tool. Run from a terminal on the VM.
set -euo pipefail
cd "$(dirname "$0")/../.."
t0=$SECONDS
step() { echo "== $1: $((SECONDS - t0)) s"; t0=$SECONDS; }

# install only what's missing, never upgrade: the Workbench's Jupyter env pins its google-cloud stack
need=$(python3 - <<'PY'
import importlib.util
mods = {"numpy": "numpy", "pandas": "pandas", "scipy": "scipy", "sklearn": "scikit-learn",
        "matplotlib": "matplotlib", "nbclient": "nbclient", "nbformat": "nbformat", "ipykernel": "ipykernel"}
print(" ".join(p for m, p in mods.items() if importlib.util.find_spec(m) is None))
PY
)
[ -z "$need" ] || python3 -m pip install -q $need

step "python"
# R packages go to R_LIBS_USER in ~ (persistent disk): slow once, seconds after. Binaries from Posit Package
# Manager when R is the system R; conda's R has to compile, so every core is used.
Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  need <- c("bigsnpr", "bigutilsr", "data.table", "dplyr", "jsonlite", "e1071")
  need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) {
    n <- parallel::detectCores(); Sys.setenv(MAKEFLAGS = paste0("-j", n))
    rel <- if (file.exists("/etc/os-release")) grep("^VERSION_CODENAME=", readLines("/etc/os-release"), value = TRUE) else character()
    rel <- sub("VERSION_CODENAME=", "", rel)
    repo <- if (length(rel) && nzchar(rel) && !grepl("conda", R.home()))
      sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", rel) else "https://cloud.r-project.org"
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                    paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
    message("installing ", paste(need, collapse = " "), " from ", repo, " on ", n, " cores")
    install.packages(need, lib = lib, repos = repo, Ncpus = n)
  }
  bioc <- c("GENESIS", "SNPRelate", "GWASTools")   # pcrelate
  bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
    BiocManager::install(bioc, lib = lib, update = FALSE, ask = FALSE, Ncpus = parallel::detectCores())
  }'
step "R packages"

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
step "plink"

make -C src
step "grm_class_tool"
command -v zstdcat >/dev/null || echo "warning: zstd missing (hm3_overlap needs it): sudo apt-get install -y zstd"
command -v gcloud >/dev/null || echo "warning: gcloud missing (aoucov sync needs it)"
df -h ~ | tail -1
echo "setup done"
