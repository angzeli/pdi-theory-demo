from __future__ import annotations
import argparse
from pathlib import Path
from .common import SYSTEMS, atom_metadata, glob_one, parse_float, read_text, read_xyz, repo_root, write_csv

def matrix(path,elements):
    started=False; cols=[]; values={}
    for line in read_text(path).splitlines():
        s=line.strip()
        if 'Bond order matrix' in s: started=True; cols=[]; continue
        if not started: continue
        f=s.split()
        if not f: continue
        if all(x.isdigit() for x in f): cols=[int(x) for x in f]; continue
        if not cols or not f[0].isdigit() or len(f[1:])!=len(cols): continue
        i=int(f[0])
        for j,raw in zip(cols,f[1:],strict=True):
            if i<j: values[(i,j)]=parse_float(raw)
    if not values: raise ValueError(f'No Mayer matrix parsed from {path}')
    return [{'atom_i':i,'element_i':elements[i],'atom_j':j,'element_j':elements[j],'bond_type':f'{elements[i]}-{elements[j]}','mayer_bond_order':v} for (i,j),v in sorted(values.items())]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); ap.add_argument('--minimum-absolute-bond-order',type=float,default=0.0); args=ap.parse_args(); repo=repo_root(args.repo); meta=atom_metadata(repo); rows=[]
    for system,cfg in SYSTEMS.items():
        path=glob_one(repo/cfg['root']/ 'bond_orders','*_mayer_bond_orders.txt'); elements={a['atom']:a['element'] for a in read_xyz(repo/cfg['xyz'])}
        for row in matrix(path,elements):
            if abs(row['mayer_bond_order'])<args.minimum_absolute_bond_order: continue
            mi=meta[system].get(row['atom_i'],{}); mj=meta[system].get(row['atom_j'],{})
            rows.append({'system':system,**row,'region_i':mi.get('region',''),'label_i':mi.get('label',''),'region_j':mj.get('region',''),'label_j':mj.get('label',''),'source_file':str(path.relative_to(repo))})
    out=repo/'results/ground_state/bond_orders/mayer_bond_orders.csv'; write_csv(out,rows); print(f'Wrote {out}'); return 0
if __name__=='__main__': raise SystemExit(main())
