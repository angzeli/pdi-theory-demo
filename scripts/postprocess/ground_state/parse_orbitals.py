from __future__ import annotations
import argparse
from pathlib import Path
from .common import SYSTEMS, parse_float, repo_root, write_csv

EV_PER_HARTREE = 27.211386245988

def molden_file(root: Path) -> Path:
    matches = sorted(p for p in root.glob("*_sp.molden.input") if p.is_file())
    if not matches:
        matches = sorted(p for p in root.glob("*.molden.input") if p.is_file())
    if not matches:
        raise FileNotFoundError(f"No Molden input file found in {root}")
    if len(matches) > 1:
        raise ValueError(f"Multiple Molden input files in {root}: {matches}")
    return matches[0]

def parse_molden(path: Path):
    rows=[]; current=None; in_mo=False
    def finish():
        if not current: return
        if "energy_au" not in current or "occupation" not in current:
            raise ValueError(f"Incomplete MO block in {path}: {current}")
        current["energy_ev"] = current["energy_au"] * EV_PER_HARTREE
        rows.append(current.copy())
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            s=line.strip()
            if s == "[MO]":
                in_mo=True; continue
            if not in_mo: continue
            if s.startswith("[") and s.endswith("]"):
                break
            if s.startswith("Sym="):
                finish(); current={"orbital_index":len(rows)+1,"symmetry":s.split("=",1)[1].strip(),"spin_channel":""}
            elif current and "=" in s:
                key,value=[part.strip() for part in s.split("=",1)]
                if key == "Ene":
                    current["energy_au"]=parse_float(value)
                elif key == "Occup":
                    current["occupation"]=parse_float(value)
                elif key == "Spin":
                    current["spin_channel"]=value
    finish()
    if not rows: raise ValueError(f"No MO blocks parsed from {path}")
    return rows

def frontier(mos):
    groups={}
    spins={row["spin_channel"] for row in mos if row["spin_channel"]}
    if len(spins) > 1:
        for row in mos: groups.setdefault(row["spin_channel"],[]).append(row)
    else:
        groups[""] = mos
    rows=[]
    for spin,group in groups.items():
        occupied=[row for row in group if row["occupation"] > 1e-8]
        virtual=[row for row in group if abs(row["occupation"]) <= 1e-8]
        if not occupied or not virtual:
            raise ValueError(f"Could not identify HOMO/LUMO for spin channel {spin or 'restricted'}")
        homo=max(occupied,key=lambda row: row["energy_au"])
        lumo=min(virtual,key=lambda row: row["energy_au"])
        row={"homo_index":homo["orbital_index"],"lumo_index":lumo["orbital_index"],"homo_energy_au":homo["energy_au"],"lumo_energy_au":lumo["energy_au"],"homo_energy_ev":homo["energy_ev"],"lumo_energy_ev":lumo["energy_ev"],"gap_ev":lumo["energy_ev"]-homo["energy_ev"]}
        if spin: row={"spin_channel":spin,**row}
        rows.append(row)
    return rows

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--repo", type=Path); args=ap.parse_args(); repo=repo_root(args.repo); frontier_rows=[]; orbital_rows=[]
    for system,cfg in SYSTEMS.items():
        path=molden_file(repo/cfg["root"]); source=str(path.relative_to(repo)); mos=parse_molden(path)
        for row in mos: orbital_rows.append({"system":system,"source_file":source,**row})
        for row in frontier(mos): frontier_rows.append({"system":system,**row,"source_file":source})
    out=repo/"results/ground_state/orbitals/frontier_orbital_energies.csv"; all_mo=repo/"results/ground_state/orbitals/orbital_energies.csv"
    frontier_fields=["system"]+(["spin_channel"] if any("spin_channel" in row for row in frontier_rows) else [])+["homo_index","lumo_index","homo_energy_au","lumo_energy_au","homo_energy_ev","lumo_energy_ev","gap_ev","source_file"]
    write_csv(out,frontier_rows,frontier_fields); write_csv(all_mo,orbital_rows,["system","source_file","orbital_index","symmetry","spin_channel","energy_au","energy_ev","occupation"]); print(f"Wrote {out}"); print(f"Wrote {all_mo}"); return 0
if __name__=="__main__": raise SystemExit(main())
