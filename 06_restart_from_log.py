#!/usr/bin/env python3
"""
Step 6 (as needed): rebuild inputs for conformers that failed, keeping everything else in the input
(routes, basis sets, GenECP blocks, title); only the link-1 geometry and OPT keyword change.

  l9999 (optimization did not converge): restart from the LAST geometry in the failed log.
  imaginary frequency (--displace):      take the optimized geometry, move it along the imaginary
      vibration (largest atom displacement 0.15 Å by default) and re-optimize, so the optimizer leaves the
      saddle point and finds the nearby minimum. Re-optimizing from the saddle point itself usually goes nowhere.

Run from the gaussian/ folder after 05_check_logs.py:
    python $PIPE_DIR/06_restart_from_log.py term_error/arcl1565_clust-7.log      # one or more logs
    python $PIPE_DIR/06_restart_from_log.py --all-l9999                          # every l9999 in check_summary.csv
    python $PIPE_DIR/06_restart_from_log.py --displace imag/arcl2321_conf-1.log  # imaginary frequency
    python $PIPE_DIR/06_restart_from_log.py --displace --all-imag                # every job in imag/

Writes <name>.com in the current folder (the original stays in term_error/ or imag/), ready to resubmit:
    ls <new>.com > jobs.txt; sbatch --array=1-$(wc -l < jobs.txt) $PIPE_DIR/04_submit_gaussian.sbatch

OPT keyword: OPT -> OPT=(MaxCycles=200,CalcFC) by default. Neither option changes the convergence criteria,
so the converged geometry is equivalent to the ArBr method. Use --opt to choose other options.
"""
import argparse, csv, re, sys
from pathlib import Path

SYM = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I", 34: "Se", 32: "Ge", 33: "As"}


def last_geometry(log):
    lines = log.read_text(errors="replace").splitlines()
    for header in ("Standard orientation:", "Input orientation:"):
        starts = [i for i, l in enumerate(lines) if header in l]
        if starts:
            i = starts[-1] + 5                            # skip the 4 header lines and the dashes
            geom = []
            while not lines[i].lstrip().startswith("-----"):
                p = lines[i].split()
                geom.append((SYM.get(int(p[1]), p[1]), float(p[3]), float(p[4]), float(p[5]))); i += 1
            return geom
    raise ValueError(f"no geometry found in {log}")


def displaced_geometry(log, amp, sign):
    """Optimized geometry moved along the lowest (imaginary) vibration from the frequency output."""
    lines = log.read_text(errors="replace").splitlines()
    f = next((i for i, l in enumerate(lines) if l.lstrip().startswith("Frequencies --")), None)
    if f is None: raise ValueError("no frequencies in log")
    freq = float(lines[f].split("--")[1].split()[0])
    if freq > 0: raise ValueError(f"lowest frequency is {freq} cm-1 (not imaginary)")
    for header in ("Standard orientation:", "Input orientation:"):    # geometry the frequencies refer to
        starts = [i for i, l in enumerate(lines[:f]) if header in l]
        if starts: break
    else: raise ValueError("no geometry before the frequencies")
    i, geom = starts[-1] + 5, []
    while not lines[i].lstrip().startswith("-----"):
        p = lines[i].split(); geom.append((SYM.get(int(p[1]), p[1]), float(p[3]), float(p[4]), float(p[5]))); i += 1
    a = next(i for i in range(f, len(lines)) if lines[i].split()[:2] == ["Atom", "AN"]) + 1
    vec = [tuple(map(float, lines[a + k].split()[2:5])) for k in range(len(geom))]
    big = max((x * x + y * y + z * z) ** 0.5 for x, y, z in vec)
    s = sign * amp / big
    return [(e, x + s * dx, y + s * dy, z + s * dz) for (e, x, y, z), (dx, dy, dz) in zip(geom, vec)], freq


def rebuild(com_text, geom, opt):
    links = com_text.split("--Link1--\n")
    first = links[0].splitlines(keepends=True)
    route_i = next(i for i, l in enumerate(first) if l.startswith("#"))
    first[route_i] = re.sub(r"\bOPT\b(=\([^)]*\)|=\S+)?", f"OPT={opt}", first[route_i], count=1, flags=re.I)
    cm_i = next(i for i in range(route_i + 1, len(first)) if re.match(r"^\s*-?\d+\s+\d+\s*$", first[i]))
    end = cm_i + 1
    while end < len(first) and first[end].strip(): end += 1
    old = [l.split()[0] for l in first[cm_i + 1:end]]
    if [g[0] for g in geom] != old:
        raise ValueError(f"atoms in log ({len(geom)}) do not match the input ({len(old)}) in the same order")
    first[cm_i + 1:end] = [f"{e:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n" for e, x, y, z in geom]
    links[0] = "".join(first)
    return "--Link1--\n".join(links)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="*")
    ap.add_argument("--all-l9999", action="store_true", help="restart every job listed with l9999 in check_summary.csv")
    ap.add_argument("--all-imag", action="store_true", help="every log in imag/ (use with --displace)")
    ap.add_argument("--displace", action="store_true", help="move along the imaginary vibration before re-optimizing")
    ap.add_argument("--amp", type=float, default=0.15, help="largest atom displacement in Å for --displace (default 0.15)")
    ap.add_argument("--minus", action="store_true", help="displace in the opposite direction (try if + returns to the saddle)")
    ap.add_argument("--opt", default="(MaxCycles=200,CalcFC)", help='OPT options, e.g. "(MaxCycles=300)"')
    a = ap.parse_args()

    logs = [Path(x) for x in a.logs]
    if a.all_l9999:
        with open("check_summary.csv") as fh:
            logs += [Path("term_error") / f"{r['job']}.log" for r in csv.DictReader(fh) if r["error_link"] == "l9999"]
    if a.all_imag:
        if not a.displace: sys.exit("--all-imag needs --displace")
        logs += sorted(Path("imag").glob("*.log"))
    if not logs: sys.exit("Give log files, or --all-l9999 / --displace --all-imag.")

    made = []
    for log in dict.fromkeys(logs):
        com = log.with_suffix(".com")
        if not com.exists(): print(f"skip {log.stem}: {com} not found"); continue
        try:
            if a.displace:
                geom, freq = displaced_geometry(log, a.amp, -1 if a.minus else 1)
                how = f"optimized geometry displaced {'-' if a.minus else '+'}{a.amp} Å along the {freq} cm-1 mode"
            else:
                geom, how = last_geometry(log), "last geometry"
            text = rebuild(com.read_text(), geom, a.opt)
        except (ValueError, StopIteration, IndexError) as e:
            print(f"skip {log.stem}: {e}"); continue
        Path(com.name).write_text(text); made.append(com.name)
        print(f"wrote {com.name}  ({how} from {log}, OPT={a.opt})")
    if made:
        Path("jobs_restart.txt").write_text("\n".join(made) + "\n")
        print(f"\n{len(made)} restart inputs. Submit with:\n  cp jobs_restart.txt jobs.txt\n"
              f"  sbatch --array=1-{len(made)} $PIPE_DIR/04_submit_gaussian.sbatch")


if __name__ == "__main__":
    main()
