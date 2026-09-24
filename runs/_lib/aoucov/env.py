"""VM setup: the plink2 binary, dsub, and a shell helper."""

import os
import subprocess

PLINK2_URL = "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
BIN_DIR = os.path.expanduser("~/bin")


def sh(script, check=True, quiet=False):
    """Run a bash snippet with `set -eo pipefail`, streaming output inline.

    Returns stdout. Used everywhere instead of hand-rolled subprocess calls so
    a failing plink step actually stops the cell.
    """
    r = subprocess.run(
        ["bash", "-c", "set -eo pipefail\n" + script],
        capture_output=True, text=True,
    )
    if not quiet and r.stdout:
        print(r.stdout.rstrip())
    if r.returncode != 0:
        print(r.stderr.rstrip())
        if check:
            raise RuntimeError(f"shell step failed ({r.returncode})")
    return r.stdout


def ensure_plink2():
    """Install plink2 into ~/bin and put it on PATH. Idempotent."""
    os.makedirs(BIN_DIR, exist_ok=True)
    if not os.access(f"{BIN_DIR}/plink2", os.X_OK):
        sh(f"""
cd /tmp
wget -q -O plink2.zip "{PLINK2_URL}"
unzip -o -q plink2.zip plink2 -d "{BIN_DIR}"
chmod +x "{BIN_DIR}/plink2"
""")
    if BIN_DIR not in os.environ["PATH"].split(":"):
        os.environ["PATH"] = f"{BIN_DIR}:{os.environ['PATH']}"
    return sh("plink2 --version", quiet=True).strip()


def ensure_dsub():
    """Install dsub and pin the cloud-sdk image its provider hard-codes.

    The pin matters: dsub ships a default cloud-sdk tag that is periodically
    withdrawn, and jobs then fail at image pull with no useful message.
    """
    sh("pip install --quiet --upgrade 'dsub>=0.5.3'")
    sh("""
DSUB_DIR=$(python -c "import dsub, os; print(os.path.dirname(dsub.__file__))")
sed -i -E "s|cloud-sdk:[0-9]+\\.[0-9]+\\.[0-9]+-slim|cloud-sdk:581.0.0-slim|g" \
  "${DSUB_DIR}/providers/google_utils.py"
dsub --version
""")
