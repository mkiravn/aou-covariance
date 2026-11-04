"""Everything here runs on a temp directory in under a second: the whole point
is that you can change a stage on a laptop and know it still works before
spending a Batch job on it."""

import json
import os
import tempfile
import unittest

from aoucov.config import RunConfig
from aoucov.manifest import pending_items, run_status, stage_status, targets
from aoucov.provenance import meta_uri, read_meta, record
from aoucov.runner import run_stage, write_task_tsv
from aoucov.store import LocalStore


def make_cfg(root):
    return RunConfig.from_dict({
        "run_id": "test",
        "root": root,
        "params": {"n_shards": 3, "seed": 1},
        "lists": {"phenotypes": ["height", "bmi"]},
        "stages": {
            "shards": {
                "dir": "03_grm", "pattern": "shard_{i}_of_{n_i}.bin",
                "over": {"i": "$n_shards"},
            },
            "cross": {
                "dir": "06_cross/{pheno}", "pattern": "acc.{i}.tsv",
                "over": {"pheno": "@phenotypes", "i": "$n_shards"},
                "needs": ["shards"],
            },
            "merged": {
                "dir": "06_cross/{pheno}", "pattern": "merged.tsv",
                "over": {"pheno": "@phenotypes"}, "needs": ["cross"],
            },
        },
    })


def touch(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


class TargetsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cfg = make_cfg(self.tmp)

    def test_counts_and_naming(self):
        t = targets(self.cfg, "shards")
        self.assertEqual(len(t), 3)
        self.assertTrue(t[0].endswith("03_grm/shard_1_of_3.bin"))

    def test_product_of_loops_is_stable(self):
        t = targets(self.cfg, "cross")
        self.assertEqual(len(t), 6)
        self.assertEqual(t, targets(self.cfg, "cross"))
        self.assertTrue(t[0].endswith("06_cross/height/acc.1.tsv"))
        self.assertTrue(t[-1].endswith("06_cross/bmi/acc.3.tsv"))

    def test_unknown_stage_names_the_options(self):
        with self.assertRaises(KeyError) as ctx:
            targets(self.cfg, "nope")
        self.assertIn("cross", str(ctx.exception))

    def test_bad_pattern_field_is_reported(self):
        cfg = make_cfg(self.tmp)
        cfg.stages["shards"].pattern = "{nosuchfield}.bin"
        with self.assertRaises(KeyError):
            targets(cfg, "shards")

    def test_config_hash_tracks_params_not_location(self):
        a = make_cfg(self.tmp)
        b = make_cfg(self.tmp)
        b.source_path = "/somewhere/else.json"
        self.assertEqual(a.config_hash(), b.config_hash())
        b.params["seed"] = 2
        self.assertNotEqual(a.config_hash(), b.config_hash())


class StatusTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cfg = make_cfg(self.tmp)
        self.store = LocalStore()

    def test_empty_run(self):
        st = stage_status(self.cfg, "shards", store=self.store)
        self.assertEqual((st.present, st.expected, st.n_missing), (0, 3, 3))
        self.assertFalse(st.done)

    def test_partial_then_complete(self):
        all_t = targets(self.cfg, "shards")
        touch(all_t[0])
        st = stage_status(self.cfg, "shards", store=self.store)
        self.assertEqual(st.n_missing, 2)
        self.assertEqual(st.missing, all_t[1:])
        for t in all_t[1:]:
            touch(t)
        self.assertTrue(stage_status(self.cfg, "shards", store=self.store).done)

    def test_empty_file_counts_as_missing(self):
        t = targets(self.cfg, "shards")[0]
        touch(t, "")
        self.assertFalse(stage_status(self.cfg, "shards", store=self.store).done)

    def test_needs_are_reported_as_blocking(self):
        for t in targets(self.cfg, "cross"):
            touch(t)
        rows = {r.stage: r for r in run_status(self.cfg, store=self.store)}
        self.assertEqual(rows["cross"].needs_unmet, ["shards"])
        self.assertEqual(rows["merged"].needs_unmet, [])

    def test_outputs_without_sidecars_are_counted(self):
        for t in targets(self.cfg, "shards"):
            touch(t)
        st = stage_status(self.cfg, "shards", store=self.store)
        self.assertEqual(st.without_meta, 3)


class ProvenanceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cfg = make_cfg(self.tmp)
        self.store = LocalStore()

    def test_sidecar_records_run_and_item(self):
        uri = targets(self.cfg, "shards")[0]
        with record(self.cfg, "shards", uri, store=self.store, item={"i": 1}) as meta:
            touch(uri)
            meta["n_pairs"] = 7
        got = read_meta(self.store, uri)
        self.assertEqual(got["run_id"], "test")
        self.assertEqual(got["stage"], "shards")
        self.assertEqual(got["item"], {"i": 1})
        self.assertEqual(got["n_pairs"], 7)
        self.assertEqual(got["config_hash"], self.cfg.config_hash())
        self.assertIn("seconds", got)
        self.assertIn("machine", got)

    def test_failure_writes_no_sidecar(self):
        uri = targets(self.cfg, "shards")[0]
        with self.assertRaises(ValueError):
            with record(self.cfg, "shards", uri, store=self.store):
                raise ValueError("job blew up")
        self.assertFalse(os.path.exists(meta_uri(uri)))

    def test_sha_from_env_is_used(self):
        os.environ["AOUCOV_GIT_SHA"] = "deadbeef"
        try:
            from aoucov import provenance
            provenance.code_version.cache_clear()
            uri = targets(self.cfg, "shards")[1]
            with record(self.cfg, "shards", uri, store=self.store):
                touch(uri)
            got = read_meta(self.store, uri)
            self.assertEqual(got["code"]["git_sha"], "deadbeef")
            self.assertEqual(got["code"]["git_sha_source"], "env")
        finally:
            del os.environ["AOUCOV_GIT_SHA"]
            from aoucov import provenance
            provenance.code_version.cache_clear()

    def test_mixed_code_is_detected(self):
        from aoucov import provenance
        for i, sha in enumerate(["aaa", "bbb", "aaa"]):
            os.environ["AOUCOV_GIT_SHA"] = sha
            provenance.code_version.cache_clear()
            uri = targets(self.cfg, "shards")[i]
            with record(self.cfg, "shards", uri, store=self.store):
                touch(uri)
        del os.environ["AOUCOV_GIT_SHA"]
        provenance.code_version.cache_clear()
        st = stage_status(self.cfg, "shards", store=self.store, check_shas=True)
        self.assertTrue(st.mixed_code)
        self.assertEqual(st.shas, {"aaa": 2, "bbb": 1})


class RunnerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cfg = make_cfg(self.tmp)
        self.store = LocalStore()
        self.calls = []

    def job(self, cfg, item, uri, meta):
        self.calls.append(item["i"])
        touch(uri, f"shard {item['i']}")
        meta["rows"] = item["i"] * 10

    def test_runs_all_then_skips(self):
        run_stage(self.cfg, "shards", self.job, store=self.store, log=lambda *a: None)
        self.assertEqual(sorted(self.calls), [1, 2, 3])
        self.calls.clear()
        run_stage(self.cfg, "shards", self.job, store=self.store, log=lambda *a: None)
        self.assertEqual(self.calls, [], "existing outputs must be skipped")

    def test_resumes_only_what_is_missing(self):
        all_t = targets(self.cfg, "shards")
        touch(all_t[0])
        touch(all_t[2])
        run_stage(self.cfg, "shards", self.job, store=self.store, log=lambda *a: None)
        self.assertEqual(self.calls, [2])

    def test_force_redoes_everything(self):
        run_stage(self.cfg, "shards", self.job, store=self.store, log=lambda *a: None)
        self.calls.clear()
        run_stage(self.cfg, "shards", self.job, store=self.store, force=True,
                  log=lambda *a: None)
        self.assertEqual(sorted(self.calls), [1, 2, 3])

    def test_parallel_matches_serial(self):
        run_stage(self.cfg, "cross", lambda c, it, u, m: touch(u, "y"),
                  jobs=4, store=self.store, log=lambda *a: None)
        self.assertTrue(stage_status(self.cfg, "cross", store=self.store).done)

    def test_dry_run_writes_nothing(self):
        run_stage(self.cfg, "shards", self.job, store=self.store, dry_run=True,
                  log=lambda *a: None)
        self.assertEqual(self.calls, [])
        self.assertEqual(stage_status(self.cfg, "shards", store=self.store).present, 0)

    def test_sidecar_written_by_runner(self):
        run_stage(self.cfg, "shards", self.job, store=self.store, log=lambda *a: None)
        st = stage_status(self.cfg, "shards", store=self.store)
        self.assertEqual(st.without_meta, 0)
        self.assertEqual(read_meta(self.store, targets(self.cfg, "shards")[2])["rows"], 30)


class TaskTsvTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cfg = make_cfg(self.tmp)
        self.tsv = os.path.join(self.tmp, "tasks.tsv")

    def test_covers_missing_only(self):
        touch(targets(self.cfg, "shards")[0])
        n = write_task_tsv(self.cfg, "shards", self.tsv)
        self.assertEqual(n, 2)
        with open(self.tsv) as fh:
            lines = fh.read().strip().split("\n")
        self.assertEqual(lines[0].split("\t"),
                         ["--env I", "--env AOUCOV_RUN", "--env AOUCOV_STAGE", "--output OUT"])
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[1].startswith("2\ttest\tshards\t"))

    def test_complete_stage_yields_no_tasks(self):
        for t in targets(self.cfg, "shards"):
            touch(t)
        self.assertEqual(write_task_tsv(self.cfg, "shards", self.tsv), 0)

    def test_multi_loop_columns(self):
        write_task_tsv(self.cfg, "cross", self.tsv)
        with open(self.tsv) as fh:
            header = fh.readline().strip().split("\t")
        self.assertEqual(header[:2], ["--env PHENO", "--env I"])


class CliTest(unittest.TestCase):
    def test_status_runs_on_a_real_config(self):
        from aoucov.__main__ import main
        tmp = tempfile.mkdtemp()
        cfg_path = os.path.join(tmp, "t.json")
        with open(cfg_path, "w") as fh:
            json.dump({
                "run_id": "t", "root": tmp,
                "params": {"n_shards": 2},
                "stages": {"s": {"dir": "d", "pattern": "f{i}.txt",
                                 "over": {"i": "$n_shards"}}},
            }, fh)
        self.assertEqual(main(["--run", cfg_path, "status"]), 0)
        self.assertEqual(main(["--run", cfg_path, "missing", "s"]), 0)


if __name__ == "__main__":
    unittest.main()
