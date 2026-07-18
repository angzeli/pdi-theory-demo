from __future__ import annotations
import argparse, re
from pathlib import Path
from .common import SYSTEMS, glob_one, read_text, repo_root, write_csv

PATTERNS = {
    "formula": [r"Formula\s*:\s*(.+)", r"Molecular formula\s*:\s*(.+)"],
    "atoms": [r"(?:Total\s+)?atoms?\s*[:=]\s*(\d+)", r"Number of atoms\s*[:=]\s*(\d+)"],
    "electrons": [r"(?:Total\s+)?electrons?\s*[:=]\s*(\d+)"],
    "basis_functions": [r"Basis functions?\s*[:=]\s*(\d+)", r"Number of basis functions\s*[:=]\s*(\d+)"],
    "occupied_orbitals": [r"Occupied orbitals?\s*[:=]\s*(\d+)"],
    "homo_index": [r"HOMO(?:\s+index)?\s*[:=]\s*(\d+)"],
    "lumo_index": [r"LUMO(?:\s+index)?\s*[:=]\s*(\d+)"],
    "multiwfn_version": [r"Multiwfn\s+([0-9][^\s,;]*)"],
}

def first(text: str, pats: list[str]) -> str:
    for p in pats:
        m = re.search(p, text, re.I)
        if m:
            return m.group(1).strip()
    return ""

def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--repo", type=Path); args = ap.parse_args()
    repo = repo_root(args.repo); rows = []
    for system, cfg in SYSTEMS.items():
        path = glob_one(repo / cfg["root"] / "wavefunction_validation", "*_validation.session.log")
        text = read_text(path); row = {"system": system, "source_file": str(path.relative_to(repo))}
        for key, pats in PATTERNS.items():
            value = first(text, pats)
            row[key] = int(value) if value and key in {"atoms","electrons","basis_functions","occupied_orbitals","homo_index","lumo_index"} else value
        low = text.lower()
        row["wavefunction_type"] = "restricted closed-shell" if "restricted closed-shell" in low else "unrestricted" if "unrestricted" in low else "restricted open-shell" if "restricted open-shell" in low else ""
        rows.append(row)
    out = repo / "results/ground_state/validation/wavefunction_summary.csv"; write_csv(out, rows); print(f"Wrote {out}"); return 0
if __name__ == "__main__": raise SystemExit(main())
