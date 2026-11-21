#!/usr/bin/env python3
"""
Step 3: write one Gaussian input (.com) per conformer, with all linked jobs.
Replaces check_manual_clustering.bash + g16_inpgen_arbrlib.py (same route lines, see method_config.py).

Run from the same folder as steps 1-2:
    python $PIPE_DIR/03_make_gaussian_inputs.py              # writes gaussian/<name>_conf-1.com, ...

File names follow the ArBr workflow (Open Babel's -m numbering): <name>_conf-1.com, <name>_conf-2.com,
or <name>_clust-1.com ... for clustered molecules.

Checks before writing (each problem goes to gaussian/problems.csv and the molecule is skipped):
  - same elements, atom count and total charge as the starting .sdf (catches a conversion that
    changed a charge or lost hydrogens)
  - even number of electrons (closed shell, multiplicity 1)
  - clustered molecules that still have a single representative are skipped (step 2 should have fixed these; see cluster_failed.txt)
"""
import argparse, csv, sys
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import method_config as M

PT = Chem.GetPeriodicTable()


def read_sdf(path):
    sup = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    mols = [m for m in sup]
    if any(m is None for m in mols): raise ValueError(f"unreadable structure in {path}")
    return mols


def summary(mol):
    return Counter(a.GetSymbol() for a in mol.GetAtoms()), sum(a.GetFormalCharge() for a in mol.GetAtoms())


def basis_block(elements, step):
    """Basis-set keyword for the route and the text block that follows the geometry (empty if not needed)."""
    heavy = sorted({e for e in elements if PT.GetAtomicNumber(e) > 36})
    light = sorted({e for e in elements if PT.GetAtomicNumber(e) <= 36}, key=PT.GetAtomicNumber)
    if not heavy:
        return (M.BASIS_OPT if step == "opt" else M.BASIS_SP), ""
    if step == "opt":
        lines = [" ".join(light) + " 0", M.BASIS_OPT, "****", " ".join(heavy) + " 0", M.HEAVY["opt_basis"], "****"]
    else:
        lines = [" ".join(light + heavy) + " 0", M.HEAVY["sp_basis"], "****"]
    lines += ["", " ".join(heavy) + " 0", M.HEAVY["ecp"]]
    return "genecp", "\n".join(lines) + "\n\n"


def make_com(stem, mol, charge, mult):
    elements = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = mol.GetConformer().GetPositions()
    geom = "".join(f"{e:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n" for e, (x, y, z) in zip(elements, xyz))
    links = []
    for L in M.LINKS:
        b_opt, block_opt = basis_block(elements, "opt")
        b_sp, block_sp = basis_block(elements, "sp")
        route = L["route"].format(basis_opt=b_opt, basis_sp=b_sp)
        block = block_opt if "{basis_opt}" in L["route"] else block_sp if "{basis_sp}" in L["route"] else ""
        text = "".join(l + "\n" for l in M.LINK0) + f"%chk={stem}.chk\n"
        text += f"{route}\n\n{stem}\n\n{charge} {mult}\n"
        text += ("" if L["read_geometry"] else geom) + "\n"
        text += block
        if L.get("nbo"): text += L["nbo"].strip() + "\n\n"
        links.append(text)
    return "--Link1--\n".join(links)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-draft", action="store_true", help="write inputs even though method_config.DRAFT is True")
    ap.add_argument("--include-single-cluster", action="store_true", help="also write molecules clustered to one representative")
    ap.add_argument("--out", default="gaussian")
    a = ap.parse_args()

    if M.DRAFT and not a.allow_draft:
        sys.exit("method_config.py is still a DRAFT. Copy the exact ArBr route lines into it and set DRAFT = False.\n"
                 "(For a single test job you can run with --allow-draft.)")
    if M.DRAFT:
        print("WARNING: using DRAFT route lines from method_config.py. Do not run production jobs with these.\n")

    out = Path(a.out); out.mkdir(exist_ok=True)
    files = sorted(Path("output").glob("*_conf.sdf")) + sorted(Path("output").glob("*_clust.sdf"))
    if not files: sys.exit("No output/*_conf.sdf or output/*_clust.sdf found. Run steps 1-2 first.")

    rows, problems, n_mol = [], [], 0
    for f in files:
        stem = f.stem                                             # arcl3949_conf or arcl3949_clust
        name = stem.rsplit("_", 1)[0]
        try:
            confs = read_sdf(f)
        except ValueError as e:
            problems.append((name, str(e))); continue
        if stem.endswith("_clust") and len(confs) == 1 and not a.include_single_cluster:
            problems.append((name, "clustered to one representative: check cluster_failed.txt (or --include-single-cluster to keep one)")); continue
        start = Path(f"{name}.sdf")
        if not start.exists():
            problems.append((name, f"starting structure {start} not found for the formula/charge check")); continue
        ref_counts, ref_charge = summary(read_sdf(start)[0])
        bad = None
        for k, m in enumerate(confs, 1):
            counts, charge = summary(m)
            if counts != ref_counts: bad = f"conformer {k}: atoms {dict(counts)} differ from start {dict(ref_counts)}"; break
            if charge != ref_charge: bad = f"conformer {k}: total charge {charge} differs from start {ref_charge}"; break
            if (sum(PT.GetAtomicNumber(e) * n for e, n in counts.items()) - charge) % 2:
                bad = "odd number of electrons (radical): not closed shell"; break
        if bad: problems.append((name, bad)); continue
        for k, m in enumerate(confs, 1):
            com = out / f"{stem}-{k}.com"
            com.write_text(make_com(com.stem, m, ref_charge, 1))
            rows.append((com.name, name, k, len(confs), ref_charge, 1, m.GetNumAtoms(),
                         int(any(PT.GetAtomicNumber(x.GetSymbol()) > 36 for x in m.GetAtoms()))))
        n_mol += 1

    with open(out / "inputs.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["com", "molecule", "conformer", "n_conformers", "charge", "multiplicity", "n_atoms", "genecp"]); w.writerows(rows)
    with open(out / "problems.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["molecule", "problem"]); w.writerows(problems)
    print(f"Wrote {len(rows)} .com files for {n_mol} molecules to {out}/ (list: {out}/inputs.csv)")
    n_gen = sum(r[-1] for r in rows)
    if n_gen:
        print(f"NOTE: {n_gen} inputs contain I and use GenECP ({M.BASIS_OPT} + {M.HEAVY['opt_basis']}/{M.HEAVY['ecp']} ECP for opt, "
              f"{M.HEAVY['sp_basis']}/{M.HEAVY['ecp']} ECP for single points). The ArBr library has no precedent for this; "
              "confirm the choice before running.")
    if problems:
        print(f"{len(problems)} molecules skipped, see {out}/problems.csv:")
        for p in problems[:20]: print("  ", *p)


if __name__ == "__main__":
    main()
