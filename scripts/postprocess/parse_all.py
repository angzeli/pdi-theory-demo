from __future__ import annotations
import argparse,subprocess,sys
from pathlib import Path
from .common import repo_root
MODULES=['build_manifest','parse_validation','parse_orbitals','parse_hirshfeld','parse_cube','parse_esp','parse_mayer','parse_qtaim','parse_dos']
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); ap.add_argument('--minimum-absolute-bond-order',type=float,default=0.0); args=ap.parse_args(); repo=repo_root(args.repo)
    for name in MODULES:
        cmd=[sys.executable,'-m',f'scripts.postprocess.{name}','--repo',str(repo)]
        if name=='parse_mayer': cmd += ['--minimum-absolute-bond-order',str(args.minimum_absolute_bond_order)]
        print('$',' '.join(cmd),flush=True); subprocess.run(cmd,cwd=repo,check=True)
    print('All parsers completed successfully.'); return 0
if __name__=='__main__': raise SystemExit(main())
