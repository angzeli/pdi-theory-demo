from __future__ import annotations
import argparse
from pathlib import Path
from .common import SYSTEMS, atom_metadata, glob_one, parse_float, read_text, repo_root, write_csv

def parse(path: Path):
    rows = []
    for line in read_text(path).splitlines():
        f = line.split()
        if not f: continue
        try:
            if len(f) >= 5 and f[0][0].isalpha():
                rows.append({"atom": len(rows)+1, "element": f[0], "x": parse_float(f[1]), "y": parse_float(f[2]), "z": parse_float(f[3]), "hirshfeld_charge": parse_float(f[4])})
            elif len(f) >= 6 and f[0].isdigit():
                rows.append({"atom": int(f[0]), "element": f[1], "x": parse_float(f[2]), "y": parse_float(f[3]), "z": parse_float(f[4]), "hirshfeld_charge": parse_float(f[5])})
        except ValueError:
            continue
    if not rows: raise ValueError(f"No charge rows parsed from {path}")
    return rows

def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--repo", type=Path); args = ap.parse_args()
    repo = repo_root(args.repo); meta = atom_metadata(repo); combined = []
    for system, cfg in SYSTEMS.items():
        path = glob_one(repo / cfg["root"] / "charges", "*_hirshfeld.chg")
        for row in parse(path):
            m = meta[system].get(int(row["atom"]), {})
            combined.append({"system": system, **row, "region": m.get("region", ""), "fragment": m.get("fragment", ""), "subregion": m.get("subregion", ""), "label": m.get("label", ""), "mapping_status": m.get("mapping_status", "unmapped"), "source_file": str(path.relative_to(repo))})
    out = repo / "results/ground_state/charges/hirshfeld_atomic_charges.csv"; write_csv(out, combined); print(f"Wrote {out}"); return 0
if __name__ == "__main__": raise SystemExit(main())
