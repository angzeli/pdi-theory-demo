#!/usr/bin/env python3
"""Parse Multiwfn Mayer bond-order matrices into a combined CSV inventory.

This script supports Multiwfn block-form bond-order matrices such as:

        1             2             3             4             5
1    3.74982331    1.03579341   -0.02186890    0.01217606   -0.01604577
2    1.03579341    3.74834249    1.25643205   -0.01593820    0.01206347

Multiwfn may print the matrix in several column blocks. This parser detects
each block automatically and reconstructs the full matrix.

Optional XYZ files can be supplied to assign element symbols to atom indices.

Example:

    python scripts/postprocess/list_mayer_bonds.py \
      --input pdi=calculations/pdi/multiwfn_analysis/bond_orders/pdi_mayer_bond_orders.txt \
      --input pdi_terminal_functionalized=calculations/pdi_terminal_functionalized/multiwfn_analysis/bond_orders/pdi_terminal_functionalized_mayer_bond_orders.txt \
      --xyz pdi=calculations/pdi/geometry_optimization/pdi_opt.xyz \
      --xyz pdi_terminal_functionalized=calculations/pdi_terminal_functionalized/geometry_optimization/pdi_terminal_functionalized_opt.xyz \
      --output results/ground_state/bond_orders/bond_inventory.csv \
      --minimum-absolute-bond-order 0.05
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True)
class MayerBond:
    """One off-diagonal Mayer bond-order matrix entry."""

    atom_i: int
    element_i: str
    atom_j: int
    element_j: str
    mayer_bond_order: float


def parse_system_path(value: str) -> tuple[str, Path]:
    """Parse a command-line value of the form SYSTEM=PATH."""

    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"Invalid value {value!r}. Expected SYSTEM=PATH."
        )

    system, raw_path = value.split("=", maxsplit=1)
    system = system.strip()
    raw_path = raw_path.strip()

    if not system:
        raise argparse.ArgumentTypeError(
            f"Invalid value {value!r}: system name is empty."
        )

    if not raw_path:
        raise argparse.ArgumentTypeError(
            f"Invalid value {value!r}: path is empty."
        )

    return system, Path(raw_path).expanduser()


def parse_float(value: str) -> float:
    """Parse ordinary or Fortran-style floating-point notation."""

    return float(value.replace("D", "E").replace("d", "e"))


def normalize_element(element: str) -> str:
    """Normalize an element symbol."""

    element = element.strip()

    if not element:
        raise ValueError("Encountered an empty element symbol.")

    return element[0].upper() + element[1:].lower()


def parse_xyz_elements(path: Path) -> dict[int, str]:
    """Return a 1-based atom-index to element-symbol mapping from an XYZ file."""

    if not path.exists():
        raise FileNotFoundError(f"XYZ file not found: {path}")

    if not path.is_file():
        raise ValueError(f"Expected an XYZ file but received: {path}")

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if not lines:
        raise ValueError(f"XYZ file is empty: {path}")

    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(
            f"Invalid XYZ atom-count line in {path}: {lines[0]!r}"
        ) from exc

    coordinate_lines = lines[2 : 2 + atom_count]

    if len(coordinate_lines) != atom_count:
        raise ValueError(
            f"XYZ file {path} declares {atom_count} atoms but contains "
            f"only {len(coordinate_lines)} coordinate lines."
        )

    elements: dict[int, str] = {}

    for atom_index, line in enumerate(coordinate_lines, start=1):
        fields = line.split()

        if len(fields) < 4:
            raise ValueError(
                f"Invalid coordinate line for atom {atom_index} "
                f"in {path}: {line!r}"
            )

        elements[atom_index] = normalize_element(fields[0])

    return elements


def is_integer_header(fields: list[str]) -> bool:
    """Return True when a line is a Multiwfn matrix column header."""

    if not fields:
        return False

    return all(field.isdigit() for field in fields)


def parse_mayer_matrix(
    path: Path,
    elements: dict[int, str] | None = None,
) -> list[MayerBond]:
    """Parse a Multiwfn block-form Mayer bond-order matrix."""

    if not path.exists():
        raise FileNotFoundError(f"Mayer bond-order file not found: {path}")

    if not path.is_file():
        raise ValueError(f"Expected a file but received: {path}")

    if path.stat().st_size == 0:
        raise ValueError(f"Mayer bond-order file is empty: {path}")

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    current_columns: list[int] = []
    matrix_started = False

    values_by_pair: dict[tuple[int, int], float] = {}
    highest_atom_index = 0

    for line in lines:
        stripped = line.strip()

        if "Bond order matrix" in stripped:
            matrix_started = True
            current_columns = []
            continue

        if not matrix_started:
            continue

        fields = stripped.split()

        if not fields:
            continue

        if is_integer_header(fields):
            current_columns = [int(field) for field in fields]
            highest_atom_index = max(
                highest_atom_index,
                max(current_columns),
            )
            continue

        if not current_columns:
            continue

        if not fields[0].isdigit():
            continue

        row_index = int(fields[0])
        raw_values = fields[1:]

        # A valid matrix row contains exactly one value per current column.
        if len(raw_values) != len(current_columns):
            continue

        highest_atom_index = max(highest_atom_index, row_index)

        for column_index, raw_value in zip(
            current_columns,
            raw_values,
            strict=True,
        ):
            try:
                bond_order = parse_float(raw_value)
            except ValueError:
                continue

            # Ignore diagonal entries because Multiwfn states that these are
            # sums of the corresponding row elements, not atom-pair bond orders.
            if row_index == column_index:
                continue

            atom_i = min(row_index, column_index)
            atom_j = max(row_index, column_index)
            key = (atom_i, atom_j)

            previous = values_by_pair.get(key)

            if previous is not None:
                # The matrix is symmetric, so the same pair may appear twice.
                if abs(previous - bond_order) > 1.0e-6:
                    raise ValueError(
                        f"Inconsistent matrix values for atoms "
                        f"{atom_i}-{atom_j} in {path}: "
                        f"{previous:.10f} versus {bond_order:.10f}."
                    )
                continue

            values_by_pair[key] = bond_order

    if not values_by_pair:
        preview_lines = [
            line.strip()
            for line in lines
            if line.strip()
        ][:20]

        preview = "\n".join(
            f"    {line}" for line in preview_lines
        )

        raise ValueError(
            f"No Mayer bond-order matrix entries could be parsed from {path}.\n"
            "First non-empty lines were:\n"
            f"{preview}"
        )

    if elements is not None:
        missing_indices = sorted(
            {
                index
                for pair in values_by_pair
                for index in pair
                if index not in elements
            }
        )

        if missing_indices:
            raise ValueError(
                f"The XYZ mapping for {path} does not contain atom indices: "
                + ", ".join(map(str, missing_indices))
            )

        if max(elements) != highest_atom_index:
            raise ValueError(
                f"Atom-count mismatch for {path}: the matrix reaches atom "
                f"{highest_atom_index}, while the XYZ file contains "
                f"{max(elements)} atoms."
            )

    bonds: list[MayerBond] = []

    for (atom_i, atom_j), bond_order in sorted(values_by_pair.items()):
        element_i = elements[atom_i] if elements else "X"
        element_j = elements[atom_j] if elements else "X"

        bonds.append(
            MayerBond(
                atom_i=atom_i,
                element_i=element_i,
                atom_j=atom_j,
                element_j=element_j,
                mayer_bond_order=bond_order,
            )
        )

    return bonds


def filter_bonds(
    bonds: Iterable[MayerBond],
    minimum_absolute_bond_order: float,
) -> list[MayerBond]:
    """Keep entries meeting the absolute Mayer bond-order threshold."""

    return [
        bond
        for bond in bonds
        if abs(bond.mayer_bond_order)
        >= minimum_absolute_bond_order
    ]


def write_inventory(
    rows: Iterable[tuple[str, MayerBond]],
    output_path: Path,
) -> None:
    """Write the combined long-form bond inventory."""

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "system",
        "atom_i",
        "element_i",
        "atom_j",
        "element_j",
        "bond_type",
        "mayer_bond_order",
    ]

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()

        for system, bond in rows:
            writer.writerow(
                {
                    "system": system,
                    "atom_i": bond.atom_i,
                    "element_i": bond.element_i,
                    "atom_j": bond.atom_j,
                    "element_j": bond.element_j,
                    "bond_type": (
                        f"{bond.element_i}-{bond.element_j}"
                    ),
                    "mayer_bond_order": (
                        f"{bond.mayer_bond_order:.10f}"
                    ),
                }
            )


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Parse one or more Multiwfn Mayer bond-order matrices "
            "and create a combined CSV inventory."
        )
    )

    parser.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="SYSTEM=PATH",
        help=(
            "Mayer bond-order file in SYSTEM=PATH form. "
            "Repeat for multiple systems."
        ),
    )

    parser.add_argument(
        "--xyz",
        action="append",
        default=[],
        metavar="SYSTEM=PATH",
        help=(
            "Optional XYZ file used to assign element symbols. "
            "Use the same system name as the corresponding --input."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output CSV path.",
    )

    parser.add_argument(
        "--minimum-absolute-bond-order",
        type=float,
        default=0.0,
        help=(
            "Keep entries whose absolute Mayer bond order is at least "
            "this value. Default: 0.0."
        ),
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the parser."""

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.minimum_absolute_bond_order < 0:
        parser.error(
            "--minimum-absolute-bond-order cannot be negative."
        )

    parsed_inputs: list[tuple[str, Path]] = []

    for value in args.input:
        try:
            parsed_inputs.append(parse_system_path(value))
        except argparse.ArgumentTypeError as exc:
            parser.error(str(exc))

    input_systems = [
        system for system, _ in parsed_inputs
    ]

    if len(input_systems) != len(set(input_systems)):
        parser.error(
            "Each --input system name must be unique."
        )

    xyz_paths: dict[str, Path] = {}

    for value in args.xyz:
        try:
            system, path = parse_system_path(value)
        except argparse.ArgumentTypeError as exc:
            parser.error(str(exc))

        if system in xyz_paths:
            parser.error(
                f"Duplicate --xyz system name: {system}"
            )

        xyz_paths[system] = path

    unknown_xyz_systems = sorted(
        set(xyz_paths) - set(input_systems)
    )

    if unknown_xyz_systems:
        parser.error(
            "--xyz contains systems not present in --input: "
            + ", ".join(unknown_xyz_systems)
        )

    combined_rows: list[tuple[str, MayerBond]] = []

    try:
        for system, mayer_path in parsed_inputs:
            elements = None

            if system in xyz_paths:
                elements = parse_xyz_elements(
                    xyz_paths[system]
                )

            bonds = parse_mayer_matrix(
                mayer_path,
                elements=elements,
            )

            unfiltered_count = len(bonds)

            bonds = filter_bonds(
                bonds,
                minimum_absolute_bond_order=(
                    args.minimum_absolute_bond_order
                ),
            )

            combined_rows.extend(
                (system, bond)
                for bond in bonds
            )

            print(
                f"Parsed {unfiltered_count} off-diagonal matrix "
                f"entries for {system!r}; retained {len(bonds)} "
                f"with |Mayer bond order| >= "
                f"{args.minimum_absolute_bond_order:g}."
            )

        combined_rows.sort(
            key=lambda item: (
                item[0],
                item[1].atom_i,
                item[1].atom_j,
            )
        )

        write_inventory(
            combined_rows,
            args.output,
        )

    except (
        FileNotFoundError,
        OSError,
        ValueError,
    ) as exc:
        print(
            f"Error: {exc}",
            file=sys.stderr,
        )
        return 1

    print(
        f"Wrote {len(combined_rows)} rows to "
        f"{args.output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())