"""dsub / Google Batch submission.

The flag set here (private address, explicit network and pet service account,
pinned cloud-sdk image) is the configuration that actually works on this
workspace. Changing any of it usually means jobs that queue forever or die at
image pull, so treat it as load-bearing.
"""

import subprocess

from .config import QC_DISK


def vcpus(machine_type):
    """n1-highmem-16 -> 16."""
    return int(machine_type.rsplit("-", 1)[-1])


def submit(run, step, script, machine_type, mem_mb, disk=QC_DISK,
           inputs=None, outputs=None, env=None):
    """Submit one Batch job and return its job id.

    `script` is a bash string; it is written to scratch and passed as --script
    rather than inlined, so quoting stays sane.
    """
    a = run.batch_args
    path = f"{run.scratch('dsub')}/{step}.sh"
    with open(path, "w") as fh:
        fh.write(script)

    flags = [
        f"--input {k}=\"{v}\"" for k, v in (inputs or {}).items()
    ] + [
        f"--output-recursive {k}=\"{v}\"" for k, v in (outputs or {}).items()
    ] + [
        f"--env {k}={v}" for k, v in {
            "VCPUS": vcpus(machine_type), "MEM_MB": mem_mb, **(env or {})
        }.items()
    ]

    cmd = f"""
dsub --provider google-batch --project {a['project']} --regions {a['region']} \
  --logging {a['logging']} --service-account {a['service_account']} \
  --network {a['network']} --subnetwork {a['subnetwork']} --use-private-address \
  --image "{a['image']}" \
  --name "{run.job_name(step)}" \
  --machine-type {machine_type} --disk-size {disk} \
  {' '.join(flags)} \
  --script {path}
"""
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    print(r.stdout or r.stderr)
    if r.returncode != 0:
        raise RuntimeError("dsub submission failed")
    return r.stdout.strip().splitlines()[-1]


def watch(run, job_id, tail=3000):
    """Current status of a job. Re-run the cell to refresh."""
    a = run.batch_args
    out = subprocess.run(["bash", "-c",
        f"dstat --provider google-batch --project {a['project']} "
        f"--location {a['region']} --jobs {job_id} --users '*' --status '*' --full"],
        capture_output=True, text=True).stdout
    print(out[-tail:])
