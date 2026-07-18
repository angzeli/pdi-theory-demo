from __future__ import annotations
import argparse, hashlib, math
from pathlib import Path
from .common import SYSTEMS, parse_float, repo_root, write_csv

def det3(a):
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def parse(path: Path):
    f=path.open('r',encoding='utf-8',errors='replace'); c1=f.readline().rstrip(); c2=f.readline().rstrip(); head=f.readline().split(); signed_n=int(head[0]); n=abs(signed_n); origin=list(map(parse_float,head[1:4])); dims=[]; axes=[]
    for _ in range(3):
        row=f.readline().split(); dims.append(abs(int(row[0]))); axes.append(list(map(parse_float,row[1:4])))
    for _ in range(n): f.readline()
    if signed_n<0:
        tokens=f.readline().split(); norb=int(tokens[0]); got=len(tokens)-1
        while got<norb: got+=len(f.readline().split())
    count=0; total=total2=totalabs=0.0; mn=math.inf; mx=-math.inf; pos=neg=zero=0
    for line in f:
        for tok in line.split():
            v=parse_float(tok); count+=1; total+=v; total2+=v*v; totalabs+=abs(v); mn=min(mn,v); mx=max(mx,v); pos+=v>0; neg+=v<0; zero+=v==0
    f.close(); expected=dims[0]*dims[1]*dims[2]
    if count!=expected: raise ValueError(f"Cube count mismatch {path}: expected {expected}, parsed {count}")
    voxel=abs(det3(axes)); mean=total/count; sd=math.sqrt(max(0,total2/count-mean*mean))
    return {"comment_1":c1,"comment_2":c2,"atom_count":n,"origin_x":origin[0],"origin_y":origin[1],"origin_z":origin[2],"nx":dims[0],"ny":dims[1],"nz":dims[2],"voxel_count":count,"voxel_volume_bohr3":voxel,"grid_volume_bohr3":voxel*count,"minimum":mn,"maximum":mx,"mean":mean,"standard_deviation":sd,"integral":total*voxel,"absolute_integral":totalabs*voxel,"positive_voxel_fraction":pos/count,"negative_voxel_fraction":neg/count,"zero_voxel_fraction":zero/count,"sha256":sha(path)}

def kind(path):
    name=path.name.lower()
    for x in ('homo','lumo','density','esp','elf','lol'):
        if x in name:return x
    return 'unknown'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path); args=ap.parse_args(); repo=repo_root(args.repo); rows=[]
    for system,cfg in SYSTEMS.items():
        root=repo/cfg['root']; paths=sorted({*root.glob('orbitals/*.cub'),*root.glob('esp/*.cub'),*root.glob('elf_lol/*.cub')})
        for path in paths: rows.append({'system':system,'cube_type':kind(path),'source_file':str(path.relative_to(repo)),**parse(path)})
    out=repo/'results/ground_state/cubes/cube_metadata.csv'; write_csv(out,rows); print(f'Wrote {out}'); return 0
if __name__=='__main__': raise SystemExit(main())
