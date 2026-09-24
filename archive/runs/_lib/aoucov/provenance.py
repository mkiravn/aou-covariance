"""Write a provenance stub next to every output, and a results line for git.

Two separate jobs:

- `write()` drops a TSV in the bucket beside the artefact. Free-form, stamped
  with the commit the code came from.
- `result_line()` prints a row to paste into the run's RESULTS.md. Keep it to
  parameters and aggregate counts -- nothing that has not cleared egress
  review goes into the repo.
"""

import datetime
import os

from .config import commit_sha


def write(out_dir, name, fields):
    path = os.path.join(out_dir, f"{name}_provenance.txt")
    with open(path, "w") as fh:
        fh.write(f"date\t{datetime.date.today().isoformat()}\n")
        fh.write(f"commit\t{commit_sha()}\n")
        for k, v in fields.items():
            fh.write(f"{k}\t{v}\n")
    print(f"provenance -> {path}")
    return path


def result_line(step, **numbers):
    """A one-line summary to copy into RESULTS.md by hand."""
    nums = "  ".join(f"{k}={v}" for k, v in numbers.items())
    line = (f"| {datetime.date.today().isoformat()} | {step} | "
            f"`{commit_sha()}` | {nums} |")
    print("\ncopy into RESULTS.md:\n" + line)
    return line
