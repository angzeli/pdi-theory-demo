from __future__ import annotations
import argparse,re
from pathlib import Path
from .common import SYSTEMS, glob_one, numbers, read_text, repo_root, write_csv

def cp_name(sig): return {'-3':'nuclear','-1':'bond','1':'ring','+1':'ring','3':'cage','+3':'cage'}.get(sig,'unknown')
def cp_signature(code): return {'1':'-3','2':'-1','3':'1','4':'3'}.get(code)

def parse_cps(path):
    rows=[]
    for ln,line in enumerate(read_text(path).splitlines(),1):
        s=line.strip(); m=re.search(r'(?:CP|critical point)\s*#?\s*(\d+).*?\(\s*3\s*,\s*([+-]?\d+)\s*\)',s,re.I)
        if m:
            vals=numbers(s); rows.append({'cp_index':int(m.group(1)),'signature':f'(3,{m.group(2)})','cp_type':cp_name(m.group(2)),'line':ln,'x':vals[-3] if len(vals)>=3 else '', 'y':vals[-2] if len(vals)>=3 else '', 'z':vals[-1] if len(vals)>=3 else '', 'raw_line':s})
    if not rows:
        for ln,line in enumerate(read_text(path).splitlines(),1):
            s=line.strip(); sig=re.search(r'\(\s*3\s*,\s*([+-]?\d+)\s*\)',s); vals=numbers(s)
            if sig and len(vals)>=4: rows.append({'cp_index':int(vals[0]),'signature':f'(3,{sig.group(1)})','cp_type':cp_name(sig.group(1)),'line':ln,'x':vals[-3],'y':vals[-2],'z':vals[-1],'raw_line':s})
    if not rows:
        for ln,line in enumerate(read_text(path).splitlines(),1):
            s=line.strip(); vals=numbers(s)
            if len(vals)==5 and str(int(vals[0]))==s.split()[0] and (sig:=cp_signature(str(int(vals[4])))):
                rows.append({'cp_index':int(vals[0]),'signature':f'(3,{sig})','cp_type':cp_name(sig),'line':ln,'x':vals[1],'y':vals[2],'z':vals[3],'raw_line':s})
    if not rows: raise ValueError(f'No critical points parsed from {path}')
    return list({r['cp_index']:r for r in rows}.values())

def parse_props(path):
    rows=[]; current=None
    for ln,line in enumerate(read_text(path).splitlines(),1):
        s=line.strip(); h=re.search(r'(?:CP|critical point)\s*#?\s*(\d+)',s,re.I)
        if h:
            if current: rows.append(current)
            current={'cp_index':int(h.group(1)),'start_line':ln,'connected_atom_i':'','connected_atom_j':'','connected_atom_i_element':'','connected_atom_j_element':''}; sig=re.search(r'\(\s*3\s*,\s*([+-]?\d+)\s*\)',s)
            if sig: current.update({'signature':f'(3,{sig.group(1)})','cp_type':cp_name(sig.group(1))})
            continue
        if current:
            connected=re.search(r'Connected atoms:\s*(\d+)\s*\(\s*([A-Za-z]+)\s*--\s*(\d+)\s*\(\s*([A-Za-z]+)',s)
            if connected:
                current.update({
                    'connected_atom_i':int(connected.group(1)),
                    'connected_atom_i_element':connected.group(2),
                    'connected_atom_j':int(connected.group(3)),
                    'connected_atom_j_element':connected.group(4),
                })
                continue
        if current and (':' in s or '=' in s):
            label=re.split(r'[:=]',s,maxsplit=1)[0].strip().lower(); vals=numbers(re.split(r'[:=]',s,maxsplit=1)[1])
            if vals:
                key=re.sub(r'[^a-z0-9]+','_',label).strip('_')
                if key and key not in current: current[key]=vals[-1]
    if current: rows.append(current)
    if not rows: raise ValueError(f'No CP properties parsed from {path}')
    return rows

def parse_paths(path):
    rows=[]; pidx=0; point=0
    for ln,line in enumerate(read_text(path).splitlines(),1):
        s=line.strip()
        if re.search(r'\bpath\b',s,re.I): pidx+=1; point=0
        vals=numbers(s)
        if len(vals)>=3:
            point+=1; rows.append({'path_index':max(pidx,1),'point_index':point,'line':ln,'x':vals[-3],'y':vals[-2],'z':vals[-1],'raw_line':s})
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); args=ap.parse_args(); repo=repo_root(args.repo); cps=[]; props=[]; paths=[]
    for system,cfg in SYSTEMS.items():
        folder=repo/cfg['root']/ 'ground_state' / 'qtaim'
        if not folder.is_dir(): folder=repo/cfg['root']/ 'qtaim'
        cp=glob_one(folder,'*_CPs.txt'); prop=glob_one(folder,'*_CPprop.txt'); path=glob_one(folder,'*_paths.txt')
        for r in parse_cps(cp): cps.append({'system':system,'source_file':str(cp.relative_to(repo)),**r})
        for r in parse_props(prop): props.append({'system':system,'source_file':str(prop.relative_to(repo)),**r})
        for r in parse_paths(path): paths.append({'system':system,'source_file':str(path.relative_to(repo)),**r})
    write_csv(repo/'results/ground_state/qtaim/critical_points.csv',cps); write_csv(repo/'results/ground_state/qtaim/critical_point_properties.csv',props); write_csv(repo/'results/ground_state/qtaim/bond_path_points.csv',paths,['system','source_file','path_index','point_index','line','x','y','z','raw_line']); print('Wrote QTAIM outputs'); return 0
if __name__=='__main__': raise SystemExit(main())
