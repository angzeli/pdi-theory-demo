from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"

SOC_ROW = re.compile(
    rf"""
    ^\s*
    (?P<triplet>\d+)\s+
    (?P<singlet>\d+)\s+
    \(\s*(?P<x_re>{FLOAT})\s*,\s*(?P<x_im>{FLOAT})\s*\)\s+
    \(\s*(?P<y_re>{FLOAT})\s*,\s*(?P<y_im>{FLOAT})\s*\)\s+
    \(\s*(?P<z_re>{FLOAT})\s*,\s*(?P<z_im>{FLOAT})\s*\)
    """,
    re.MULTILINE | re.VERBOSE,
)


def as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def extract_soc_section(text: str) -> str:
    marker = "CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS"

    if marker not in text:
        raise ValueError("SOC matrix-element section was not found.")

    section = text.split(marker, 1)[1]

    possible_end_markers = (
        "SOC stabilization of the ground state",
        "CALCULATED ABSORPTION SPECTRUM",
        "ORCA TERMINATED NORMALLY",
    )

    end_positions = [
        section.find(end_marker)
        for end_marker in possible_end_markers
        if section.find(end_marker) >= 0
    ]

    if end_positions:
        section = section[: min(end_positions)]

    return section


def parse_soc(path: Path) -> list[dict[str, float | int]]:
    text = path.read_text(errors="replace")

    if "ORCA TERMINATED NORMALLY" not in text:
        raise ValueError(f"{path} did not terminate normally.")

    section = extract_soc_section(text)
    rows: list[dict[str, float | int]] = []

    for match in SOC_ROW.finditer(section):
        values = {
            key: as_float(match.group(key))
            for key in ("x_re", "x_im", "y_re", "y_im", "z_re", "z_im")
        }

        magnitude = math.sqrt(sum(value**2 for value in values.values()))

        rows.append(
            {
                "triplet": int(match.group("triplet")),
                "singlet": int(match.group("singlet")),
                **values,
                "soc_cm-1": magnitude,
            }
        )

    if not rows:
        raise ValueError(
            "The SOC section was found, but no SOC rows matched the parser. "
            "Inspect the ORCA table formatting."
        )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--csv", type=Path)
    parser.add_argument(
        "--singlet",
        type=int,
        default=1,
        help="Only print couplings involving this singlet root.",
    )
    args = parser.parse_args()

    rows = parse_soc(args.output)
    selected = [row for row in rows if row["singlet"] == args.singlet]
    selected.sort(key=lambda row: row["triplet"])

    print(f"\nSOC values involving S{args.singlet}")
    print("Pair       SOC / cm^-1")
    print("----------------------")

    for row in selected:
        print(
            f"S{row['singlet']}-T{row['triplet']:<3d} "
            f"{row['soc_cm-1']:12.6f}"
        )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)

        with args.csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()