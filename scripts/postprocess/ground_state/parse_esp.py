from __future__ import annotations
import argparse,re,sys
from pathlib import Path
from .common import SYSTEMS, numbers, read_text, repo_root, write_csv

def parse_session_metrics(path):
    metrics=[]
    keys=('minimum','maximum','average','variance','standard deviation','molecular polarity index','positive area','negative area','surface area','polar surface','nonpolar surface','balance','internal charge separation','skewness')
    for ln,line in enumerate(read_text(path).splitlines(),1):
        s=line.strip(); low=s.lower()
        if ':' not in s and '=' not in s: continue
        label,rhs=re.split(r'[:=]',s,maxsplit=1); vals=numbers(rhs)
        if s and any(k in low for k in keys) and vals:
            metrics.append({'line':ln,'label':label.strip(),'value':vals[0],'raw_line':s})
    return metrics

def parse_stats(path):
    extrema=[]; section=None
    for ln,line in enumerate(read_text(path).splitlines(),1):
        s=line.strip(); low=s.lower()
        if 'number of surface minima' in low:
            section='minimum'; continue
        if 'number of surface maxima' in low:
            section='maximum'; continue
        vals=numbers(s)
        if section and len(vals)>=7:
            extrema.append({'line':ln,'extremum_type':section,'x':vals[4],'y':vals[5],'z':vals[6],'esp_value':vals[1],'raw_line':s,'is_global':s.startswith('*')})
    if not extrema: raise ValueError(f'No ESP extrema parsed from {path}')
    return extrema

def derived_global_metrics(extrema):
    metrics=[]
    for kind,selector in [('minimum',min),('maximum',max)]:
        rows=[r for r in extrema if r['extremum_type']==kind]
        if not rows: continue
        marked=[r for r in rows if r.get('is_global')]
        row=marked[0] if marked else selector(rows,key=lambda r:r['esp_value'])
        metrics.append({'line':row['line'],'label':f'Derived global surface {kind} (a.u.)','value':row['esp_value'],'raw_line':row['raw_line']})
    return metrics

def parse_pdb(path):
    rows=[]
    for line in read_text(path).splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        try: rows.append({'record':line[:6].strip(),'serial':int(line[6:11]),'atom_name':line[12:16].strip(),'residue_name':line[17:20].strip(),'x':float(line[30:38]),'y':float(line[38:46]),'z':float(line[46:54]),'occupancy':float(line[54:60] or 0),'beta_factor':float(line[60:66] or 0),'element':line[76:78].strip()})
        except ValueError: continue
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); args=ap.parse_args(); repo=repo_root(args.repo); metrics=[]; text_ext=[]; pdb_ext=[]
    for system,cfg in SYSTEMS.items():
        folder=repo/cfg['root']/ 'esp'; stats=sorted(folder.glob('*_esp_surface_statistics.txt')); pdbs=sorted(folder.glob('*_esp_extrema.pdb')); logs=sorted(folder.glob('*_esp_surface.session.log'))
        if not stats or not pdbs: raise FileNotFoundError(f'Missing ESP statistics/extrema in {folder}')
        e=parse_stats(stats[0]); m=parse_session_metrics(logs[0]) if logs else []; m_source=logs[0] if m and logs else stats[0]
        if not m:
            print(f'Warning: no ESP surface summary metrics found for {system}; deriving only global extrema from {stats[0]}',file=sys.stderr)
            m=derived_global_metrics(e)
        for r in m: metrics.append({'system':system,'source_file':str(m_source.relative_to(repo)),**r})
        for r in e: text_ext.append({'system':system,'source_file':str(stats[0].relative_to(repo)),**r})
        for r in parse_pdb(pdbs[0]): pdb_ext.append({'system':system,'source_file':str(pdbs[0].relative_to(repo)),**r})
    write_csv(repo/'results/ground_state/esp/esp_surface_metrics.csv',metrics); write_csv(repo/'results/ground_state/esp/esp_extrema_from_text.csv',text_ext,['system','source_file','line','extremum_type','x','y','z','esp_value','raw_line']); write_csv(repo/'results/ground_state/esp/esp_extrema_from_pdb.csv',pdb_ext,['system','source_file','record','serial','atom_name','residue_name','x','y','z','occupancy','beta_factor','element']); print('Wrote ESP outputs'); return 0
if __name__=='__main__': raise SystemExit(main())
