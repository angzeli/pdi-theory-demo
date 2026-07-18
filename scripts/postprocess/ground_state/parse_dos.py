from __future__ import annotations
import argparse,re
from pathlib import Path
from .common import SYSTEMS, glob_one, parse_float, read_text, repo_root, write_csv, write_json

def table(path):
    rows=[]; width=None
    for line in read_text(path).splitlines():
        try: vals=[parse_float(x) for x in line.split()]
        except ValueError: continue
        if not vals: continue
        width=width or len(vals)
        if len(vals)==width: rows.append(vals)
    if not rows: raise ValueError(f'No numeric data in {path}')
    return rows

def settings(path,session):
    text=(read_text(path)+'\n' if path else '')+(read_text(session) if session else ''); out={"energy_unit":"eV" if re.search(r'\beV\b',text) else '',"broadening":"Gaussian" if re.search('Gaussian',text,re.I) else '',"projection_method":"Mulliken" if re.search('Mulliken',text,re.I) else ''}
    m=re.search(r'Energy range\s*:\s*([+-]?\d+(?:\.\d+)?)\s*(?:to|[-–])\s*([+-]?\d+(?:\.\d+)?)',text,re.I)
    out['energy_min_ev']=float(m.group(1)) if m else ''; out['energy_max_ev']=float(m.group(2)) if m else ''
    m=re.search(r'Step\s*:\s*([0-9.]+)',text,re.I); out['energy_step_ev']=float(m.group(1)) if m else ''
    vals=re.findall(r'FWHM\s*[:=]\s*([0-9.]+)\s*eV',text,re.I); out['fwhm_ev']=float(vals[-1]) if vals else ''
    return out

def convert(system,source,rows):
    width=len(rows[0]); names=['energy_ev','tdos','opdos']+[f'pdos_fragment_{i}' for i in range(1,max(1,width-2))]
    names=names[:width]; return [{'system':system,'source_file':source,**dict(zip(names,r,strict=True))} for r in rows]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); args=ap.parse_args(); repo=repo_root(args.repo); curves=[]; lines=[]; meta={}
    for system,cfg in SYSTEMS.items():
        folder=repo/cfg['root']/ 'dos'; curve=glob_one(folder,'*_dos_curve.txt'); line=glob_one(folder,'*_dos_line.txt'); frag=glob_one(folder,'*_DOSfrag.txt'); session=glob_one(folder,'*_dos.session.log',required=False); setting=glob_one(folder,'*_dos_settings.txt',required=False)
        curves+=convert(system,str(curve.relative_to(repo)),table(curve)); lines+=convert(system,str(line.relative_to(repo)),table(line)); meta[system]={**settings(setting,session),'curve_file':str(curve.relative_to(repo)),'line_file':str(line.relative_to(repo)),'fragment_file':str(frag.relative_to(repo))}
    write_csv(repo/'results/ground_state/dos/dos_curves.csv',curves); write_csv(repo/'results/ground_state/dos/dos_lines.csv',lines); write_json(repo/'results/ground_state/dos/dos_metadata.json',meta); print('Wrote DOS outputs'); return 0
if __name__=='__main__': raise SystemExit(main())
