from __future__ import annotations

import csv, json, re
from pathlib import Path
from typing import Any, Iterable

SYSTEMS = {
    "pdi": {
        "root": Path("calculations/pdi/multiwfn_analysis"),
        "xyz": Path("calculations/pdi/electronic_structure_calculation/pdi_opt.xyz"),
    },
    "pdi_terminal_functionalized": {
        "root": Path("calculations/pdi_terminal_functionalized/multiwfn_analysis"),
        "xyz": Path("calculations/pdi_terminal_functionalized/electronic_structure_calculation/pdi_terminal_functionalized_opt.xyz"),
    },
}


def repo_root(start: Path | None = None) -> Path:
    p = (start or Path.cwd()).resolve()
    if p.is_file():
        p = p.parent
    for candidate in [p, *p.parents]:
        if (candidate / "pyproject.toml").is_file():
            return candidate
    raise FileNotFoundError("Could not find repository root containing pyproject.toml")


def read_text(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")
    return path.read_text(encoding="utf-8", errors="replace")


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            raise ValueError(f"No rows and no fieldnames for {path}")
        fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def glob_one(folder: Path, pattern: str, required: bool = True) -> Path | None:
    matches = sorted(p for p in folder.glob(pattern) if p.is_file())
    if not matches:
        if required:
            raise FileNotFoundError(f"No {pattern!r} in {folder}")
        return None
    if len(matches) > 1:
        raise ValueError(f"Multiple matches for {pattern!r} in {folder}: {matches}")
    return matches[0]


def read_xyz(path: Path) -> list[dict[str, Any]]:
    lines = read_text(path).splitlines()
    n = int(lines[0].strip())
    coords = lines[2:2+n]
    if len(coords) != n:
        raise ValueError(f"XYZ count mismatch: {path}")
    rows = []
    for i, line in enumerate(coords, 1):
        f = line.split()
        rows.append({"atom": i, "element": f[0], "x": parse_float(f[1]), "y": parse_float(f[2]), "z": parse_float(f[3])})
    return rows


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def atom_metadata(repo: Path) -> dict[str, dict[int, dict[str, str]]]:
    out = {"pdi": {}, "pdi_terminal_functionalized": {}}
    mapping = repo / "config/pdi_core_atom_mapping.csv"
    if mapping.is_file():
        for row in read_csv_dicts(mapping):
            p = row.get("parent_atom", "").strip()
            f = row.get("functionalized_atom", "").strip()
            base = {"region": row.get("region", "").strip(), "label": row.get("label", "").strip()}
            if p:
                out["pdi"][int(p)] = {**base, "mapping_status": "matched" if f else "parent_only"}
            if f:
                out["pdi_terminal_functionalized"][int(f)] = {**base, "mapping_status": "matched" if p else "functionalized_only"}
    regions = repo / "config/functionalized_atom_regions.csv"
    if regions.is_file():
        for row in read_csv_dicts(regions):
            system = row.get("system", "").strip()
            if system not in out:
                continue
            atom = int(row["atom"])
            fragment = row.get("fragment", "").strip()
            subregion = row.get("subregion", "").strip()
            region = row.get("region", "").strip() or (f"{fragment}_{subregion}" if fragment and subregion else fragment)
            out[system][atom] = {
                **out[system].get(atom, {}),
                "region": region,
                "fragment": fragment,
                "subregion": subregion,
                "label": row.get("label", "").strip(),
                "mapping_status": out[system].get(atom, {}).get("mapping_status", "functionalized_only"),
            }
    return out


def numbers(text: str) -> list[float]:
    pat = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    return [parse_float(x) for x in re.findall(pat, text)]
