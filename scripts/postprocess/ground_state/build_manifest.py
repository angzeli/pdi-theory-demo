from __future__ import annotations
import argparse,hashlib
from pathlib import Path
from .common import SYSTEMS, repo_root, write_csv

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); args=ap.parse_args(); repo=repo_root(args.repo); rows=[]
    for system,cfg in SYSTEMS.items():
        root=repo/cfg['root']
        for path in sorted(root.rglob('*')):
            if path.is_file() and path.name!='.DS_Store': rows.append({'system':system,'analysis_folder':path.parent.name,'filename':path.name,'relative_path':str(path.relative_to(repo)),'size_bytes':path.stat().st_size,'sha256':digest(path)})
    out=repo/'results/ground_state/manifest/input_manifest.csv'; write_csv(out,rows); print(f'Wrote {out}'); return 0
if __name__=='__main__': raise SystemExit(main())
