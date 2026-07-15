#!/usr/bin/env python3

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path


def compress_ranges(indices: list[int]) -> str:
    if not indices:
        return ""

    ranges: list[str] = []
    start = previous = indices[0]

    for value in indices[1:]:
        if value == previous + 1:
            previous = value
            continue

        ranges.append(
            str(start) if start == previous else f"{start}-{previous}"
        )
        start = previous = value

    ranges.append(
        str(start) if start == previous else f"{start}-{previous}"
    )

    return ",".join(ranges)


def read_atoms(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = defaultdict(list)
    in_atoms = False

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if line.lower().startswith("[atoms]"):
                in_atoms = True
                continue

            if in_atoms and line.startswith("["):
                break

            if not in_atoms or not line:
                continue

            fields = line.split()

            if len(fields) < 2:
                continue

            symbol = fields[0]
            atom_index = int(fields[1])
            groups[symbol].append(atom_index)

    return dict(groups)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: list_molden_atoms.py FILE.molden.input"
        )

    path = Path(sys.argv[1]).expanduser().resolve()

    if not path.is_file():
        raise SystemExit(f"File not found: {path}")

    groups = read_atoms(path)

    print(f"File: {path}")
    print()

    for element, indices in groups.items():
        print(f"{element:>2}: {compress_ranges(indices)}")


if __name__ == "__main__":
    main()
