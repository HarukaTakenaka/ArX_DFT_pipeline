#!/usr/bin/env python3
"""
Step 5: check finished Gaussian jobs and sort them. Replaces post_processing_v1_1_0.bash
+ log_check_for_processing_v1_1_0.py, with the same sorting:

    logs/        finished normally, no imaginary frequency           (.log moved here)
    term_error/  wrong number of "Normal termination" lines            (.log and .com moved here)
    imag/        lowest frequency <= 0                                (.log and .com moved here)
    resubmit/    .com with no .log                                    (.com moved here)

Run from the gaussian/ folder AFTER all jobs have finished (a log that is still being written
looks like a termination error):
    python $PIPE_DIR/05_check_logs.py            # expected terminations taken from method_config.py
    python $PIPE_DIR/05_check_logs.py --clean    # also delete *.chk, *.fchk, fort.7, slurm-* (the old script always did)

Changes from the old scripts:
  - the expected number of "Normal termination" lines is read from method_config.py (or --expected),
    instead of being edited by hand on line 36
  - frequencies are read from the "Frequencies --" line, not fixed character columns
  - writes check_summary.csv (one row per conformer, with the Gaussian error link, e.g. l9999 = optimization
    did not converge, l502 = SCF did not converge) and molecule_status.csv (complete only if every
    conformer of the molecule is in logs/)
  - safe to rerun: new logs are checked, files already sorted are left alone
"""
import argparse, csv, re, shutil, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import method_config as M
    DEFAULT_EXPECTED = M.EXPECTED_NORMAL_TERMINATIONS
except Exception:
    DEFAULT_EXPECTED = 4


def inspect(log):
    text = log.read_text(errors="replace").splitlines()
    n_normal = sum("Normal termination" in l for l in text)
    err = next((l.strip() for l in reversed(text) if "Error termination" in l), "")
    link = (re.search(r"(l\d{3,4})\.exe", err) or [None, ""])[1] if err else ""
    freqs = next((l for l in text if l.lstrip().startswith("Frequencies --")), None)
    low = min(float(x) for x in freqs.split("--")[1].split()) if freqs else None
    return n_normal, err, link, low


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expected", type=int, default=DEFAULT_EXPECTED, help="Normal termination lines per finished log")
    ap.add_argument("--clean", action="store_true", help="delete *.chk, *.fchk, fort.7 and slurm-* files first")
    a = ap.parse_args()

    here = Path(".")
    if a.clean:
        for pat in ("*.chk", "*.fchk", "fort.7", "slurm-*"):
            for f in here.glob(pat): f.unlink()
    for d in ("logs", "term_error", "imag", "resubmit"): (here / d).mkdir(exist_ok=True)

    rows = []
    for com in sorted(here.glob("*.com")):
        log = com.with_suffix(".log")
        if not log.exists() and (here / "logs" / log.name).exists():
            continue                                   # finished and sorted on an earlier run (the .com stays here)
        if not log.exists():
            shutil.move(com, here / "resubmit" / com.name)
            rows.append((com.stem, "resubmit", "", "", "", "no log")); continue
        n, err, link, low = inspect(log)
        if n != a.expected:
            status = "term_error"
            for f in (log, com): shutil.move(f, here / status / f.name)
        elif low is None or low <= 0:
            status = "imag"
            for f in (log, com): shutil.move(f, here / status / f.name)
        else:
            status = "ok"
            shutil.move(log, here / "logs" / log.name)
        if status == "term_error" and not err:
            note = "no Gaussian error message: job was killed (preempted or hit the time limit); resubmit as is"
        else:
            note = err if status == "term_error" else ("no frequencies found" if low is None else "")
        rows.append((com.stem, status, n, link, "" if low is None else low, note))

    # previously sorted logs count towards molecule status too
    seen = {r[0] for r in rows}
    for d, status in (("logs", "ok"), ("term_error", "term_error"), ("imag", "imag")):
        for f in (here / d).glob("*.log"):
            if f.stem not in seen: rows.append((f.stem, status, "", "", "", "sorted earlier")); seen.add(f.stem)
    for f in (here / "resubmit").glob("*.com"):
        if f.stem not in seen: rows.append((f.stem, "resubmit", "", "", "", "sorted earlier")); seen.add(f.stem)

    with open("check_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["job", "status", "normal_terminations", "error_link", "lowest_freq_cm-1", "note"]); w.writerows(rows)

    mol = defaultdict(list)
    for job, status, *_ in rows: mol[re.sub(r"_(conf|clust)-\d+$", "", job)].append(status)
    with open("molecule_status.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["molecule", "n_conformers", "n_ok", "complete"])
        for m, s in sorted(mol.items()): w.writerow([m, len(s), s.count("ok"), all(x == "ok" for x in s)])

    count = lambda s: sum(r[1] == s for r in rows)
    print(f"Expected {a.expected} 'Normal termination' lines per log.")
    print(f"ok {count('ok')} | term_error {count('term_error')} | imag {count('imag')} | resubmit {count('resubmit')}")
    print(f"Molecules complete: {sum(all(x == 'ok' for x in s) for s in mol.values())}/{len(mol)} (molecule_status.csv)")
    for r in rows:
        if r[1] == "imag" and r[5] != "sorted earlier": print(f"  imaginary: {r[0]}  {r[4]} cm-1")
        if r[1] == "term_error" and r[5] != "sorted earlier":
            print(f"  error:     {r[0]}  {r[2]} terminations  {r[3] or 'killed (no error message)'}")


if __name__ == "__main__":
    main()
