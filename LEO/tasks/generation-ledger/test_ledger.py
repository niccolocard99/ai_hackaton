"""Tests for ledger.py. Run from this folder:  python3 -m unittest test_ledger -v

Everything happens in temporary folders; the real runs.csv and variants.csv are never touched.
"""
import contextlib
import csv
import io
import os
import tempfile
import unittest
from pathlib import Path

import ledger

HEAVY_A = "EVQLVESGGGLVQPGGSLRLSCAAS"
HEAVY_B = "QVQLQQSGAELARPGASVKMSCKAS"
LIGHT = "DIQMTQSPSSLSASVGDRVTITC"


def run(*argv):
    """Call ledger.main in process; return (exit_code, stdout)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        try:
            code = ledger.main([str(x) for x in argv])
        except SystemExit as e:
            code = 1 if isinstance(e.code, str) else e.code
            if isinstance(e.code, str):
                out.write(e.code)
    return code, out.getvalue()


class LedgerTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)
        self.led = self.tmp / "ledger"
        self.led.mkdir()
        self.spec = self.tmp / "spec.yaml"
        self.cif = self.tmp / "target.cif"
        self.cif.write_text("data_target\n")
        self.spec.write_text("design:\n  input: target.cif\n  conditions: {}\n")

    def runflags(self, **kw):
        flags = ["--ledger-dir", self.led, "--method", "design", "--model",
                 "quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09",
                 "--input", self.spec, "--seed", "1", "--setting", "num-replicates=100"]
        for k, v in kw.items():
            flags += [f"--{k.replace('_', '-')}", v]
        return flags

    def rows(self, name):
        with open(self.led / name, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))

    def designs(self, rows, name="designs.csv", cols=("sequence_id", "kind", "heavy_sequence",
                                                      "light_sequence", "cdrh3", "base_validation")):
        path = self.tmp / name
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(cols)
            w.writerows(rows)
        return path

    # ---- runs ----

    def test_new_run_then_seen(self):
        self.assertEqual(run("check-run", *self.runflags())[0], 0)
        self.assertEqual(run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")[0], 0)
        code, out = run("check-run", *self.runflags())
        self.assertEqual(code, 3)
        self.assertIn("r1", out)

    def test_model_spelling_is_normalised(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        flags = self.runflags()
        flags[flags.index("--model") + 1] = "apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09"
        self.assertEqual(run("check-run", *flags)[0], 3)

    def test_changed_input_content_is_a_new_run(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        self.spec.write_text("design:\n  input: target.cif\n  conditions: {x: 1}\n")
        self.assertEqual(run("check-run", *self.runflags())[0], 0)

    def test_file_the_spec_points_at_is_hashed(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        self.cif.write_text("data_target\nATOM changed\n")  # spec untouched, structure edited
        self.assertEqual(run("check-run", *self.runflags())[0], 0)

    def test_seed_and_setting_change_the_run(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        flags = self.runflags()
        flags[flags.index("--seed") + 1] = "2"
        self.assertEqual(run("check-run", *flags)[0], 0)
        flags = self.runflags()
        flags[flags.index("--setting") + 1] = "num-replicates=200"
        self.assertEqual(run("check-run", *flags)[0], 0)

    def test_setting_order_does_not_matter(self):
        a = self.runflags() + ["--setting", "sequence-tsr-k=2"]
        run("add-run", *a, "--run-id", "r1", "--owner", "t")
        b = self.runflags()
        b[b.index("--setting") + 1] = "sequence-tsr-k=2"
        b += ["--setting", "num-replicates=100"]
        self.assertEqual(run("check-run", *b)[0], 3)

    def test_output_path_is_refused_as_a_setting(self):
        code, out = run("check-run", *self.runflags(), "--setting", "output=runs/x")
        self.assertEqual(code, 1)
        self.assertIn("output path", out)

    def test_add_run_refuses_a_blocking_duplicate_and_a_reused_id(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        self.assertEqual(run("add-run", *self.runflags(), "--run-id", "r2", "--owner", "t")[0], 1)
        self.assertEqual(run("add-run", *self.runflags(seed="9"), "--run-id", "r1", "--owner", "t")[0], 1)
        self.assertEqual(len(self.rows("runs.csv")), 1)

    def test_failed_run_may_be_retried_but_empty_may_not(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        run("set-status", "--ledger-dir", self.led, "--run-id", "r1", "--status", "failed", "--notes", "OOM")
        code, out = run("check-run", *self.runflags())
        self.assertEqual(code, 0)
        self.assertIn("FAILED", out)
        run("set-status", "--ledger-dir", self.led, "--run-id", "r1", "--status", "empty")
        self.assertEqual(run("check-run", *self.runflags())[0], 3)
        self.assertEqual(len(self.rows("runs.csv")), 3)  # append-only: one row per change

    def test_missing_input_is_an_error(self):
        flags = self.runflags()
        flags[flags.index("--input") + 1] = self.tmp / "nope.yaml"
        self.assertEqual(run("check-run", *flags)[0], 1)

    # ---- variants ----

    def test_variants_new_seen_and_dup_in_file(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        first = self.designs([("d1", "vhh", HEAVY_A, "", "", ""), ("d2", "vhh", HEAVY_B, "", "", "")])
        code, out = run("add", first, "--ledger-dir", self.led, "--run-id", "r1")
        self.assertEqual(code, 0)
        self.assertEqual(len(self.rows("variants.csv")), 2)

        second = self.designs([("x1", "vhh", HEAVY_A.lower(), "", "", ""),  # case does not matter
                               ("x2", "vhh", HEAVY_A, "", "", ""),           # repeated in the file
                               ("x3", "vhh", HEAVY_A + "W", "", "", "")])    # genuinely new
        fresh = self.tmp / "fresh.csv"
        code, out = run("check", second, "--ledger-dir", self.led, "--new-only", fresh)
        self.assertEqual(code, 3)
        self.assertRegex(out, r"SEEN\s+x1")
        self.assertRegex(out, r"DUP-IN-FILE\s+x2")
        self.assertRegex(out, r"NEW\s+x3")
        with open(fresh, newline="") as f:
            self.assertEqual([r["sequence_id"] for r in csv.DictReader(f)], ["x3"])

    def test_light_chain_is_part_of_the_identity(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        run("add", self.designs([("p1", "paired", HEAVY_A, LIGHT, "", "")]), "--ledger-dir", self.led,
            "--run-id", "r1")
        code, _ = run("check", self.designs([("p2", "paired", HEAVY_A, "", "", "")], "b.csv"),
                      "--ledger-dir", self.led)
        self.assertEqual(code, 0)  # same heavy, no light: a different antibody

    def test_same_cdrh3_is_reported_but_not_blocking(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        run("add", self.designs([("d1", "vhh", HEAVY_A, "", "ARDYYGSSYFDY", "")]), "--ledger-dir", self.led,
            "--run-id", "r1")
        code, out = run("check", self.designs([("d9", "vhh", HEAVY_B, "", "ARDYYGSSYFDY", "")], "c.csv"),
                        "--ledger-dir", self.led)
        self.assertEqual(code, 0)
        self.assertIn("CDRH3 ARDYYGSSYFDY is already in the ledger", out)

    def test_add_is_idempotent_and_a_verdict_update_is_a_new_row(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        table = self.designs([("d1", "vhh", HEAVY_A, "", "", "")])
        run("add", table, "--ledger-dir", self.led, "--run-id", "r1")
        run("add", table, "--ledger-dir", self.led, "--run-id", "r1")
        self.assertEqual(len(self.rows("variants.csv")), 1)
        table = self.designs([("d1", "vhh", HEAVY_A, "", "ARDYYGSSYFDY", "pass")])
        run("add", table, "--ledger-dir", self.led, "--run-id", "r1")
        rows = self.rows("variants.csv")
        self.assertEqual([r["base_validation"] for r in rows], ["untested", "PASS"])
        code, out = run("check", table, "--ledger-dir", self.led)
        self.assertEqual(code, 3)
        self.assertIn("validation=PASS", out)

    def test_add_needs_a_registered_run_and_a_valid_verdict(self):
        table = self.designs([("d1", "vhh", HEAVY_A, "", "", "")])
        self.assertEqual(run("add", table, "--ledger-dir", self.led, "--run-id", "ghost")[0], 1)
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        bad = self.designs([("d1", "vhh", HEAVY_A, "", "", "maybe")], "bad.csv")
        self.assertEqual(run("add", bad, "--ledger-dir", self.led, "--run-id", "r1")[0], 1)
        self.assertFalse((self.led / "variants.csv").exists())

    def test_bad_sequence_is_rejected(self):
        bad = self.designs([("d1", "vhh", "EVQL*", "", "", "")], "bad.csv")
        self.assertEqual(run("check", bad, "--ledger-dir", self.led)[0], 1)

    def test_edited_header_is_refused(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        text = (self.led / "runs.csv").read_text()
        (self.led / "runs.csv").write_text(text.replace("run_key", "key", 1))
        self.assertEqual(run("check-run", *self.runflags())[0], 1)

    def test_missing_final_newline_does_not_glue_rows(self):
        run("add-run", *self.runflags(), "--run-id", "r1", "--owner", "t")
        path = self.led / "runs.csv"
        path.write_bytes(path.read_bytes().rstrip(b"\n"))
        run("add-run", *self.runflags(seed="2"), "--run-id", "r2", "--owner", "t")
        self.assertEqual([r["run_id"] for r in self.rows("runs.csv")], ["r1", "r2"])


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent)
    unittest.main()
