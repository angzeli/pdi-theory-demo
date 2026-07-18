from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


OUTPUT_COLUMNS = [
    "bond_id",
    "label",
    "parent_i",
    "parent_j",
    "functionalized_i",
    "functionalized_j",
    "comparison_type",
]


def canonical_pair(atom_i: int, atom_j: int) -> tuple[int, int]:
    """Return an atom pair in ascending numerical order."""
    return tuple(sorted((int(atom_i), int(atom_j))))


def element_pair(element_i: str, element_j: str) -> str:
    """Return a readable element-pair label."""
    return f"{element_i}-{element_j}"


def prepare_inventory(
    inventory: pd.DataFrame,
    minimum_bond_order: float,
    include_hydrogen: bool,
) -> pd.DataFrame:
    """Filter and canonicalize parsed Mayer bond-order records."""

    required = {
        "system",
        "atom_i",
        "element_i",
        "atom_j",
        "element_j",
        "mayer_bond_order",
    }

    missing = required - set(inventory.columns)
    if missing:
        raise ValueError(
            "Bond inventory is missing required columns: "
            + ", ".join(sorted(missing))
        )

    data = inventory.copy()

    data["atom_i"] = data["atom_i"].astype(int)
    data["atom_j"] = data["atom_j"].astype(int)
    data["mayer_bond_order"] = pd.to_numeric(
        data["mayer_bond_order"],
        errors="coerce",
    )

    data = data.dropna(subset=["mayer_bond_order"])

    # For a selected-bond file, retain chemically bonded positive interactions.
    data = data[
        data["mayer_bond_order"] >= minimum_bond_order
    ].copy()

    if not include_hydrogen:
        data = data[
            (data["element_i"] != "H")
            & (data["element_j"] != "H")
        ].copy()

    pairs = data.apply(
        lambda row: canonical_pair(row["atom_i"], row["atom_j"]),
        axis=1,
    )

    data["pair_i"] = [pair[0] for pair in pairs]
    data["pair_j"] = [pair[1] for pair in pairs]

    # Ensure element ordering follows atom ordering.
    swap_mask = data["atom_i"] > data["atom_j"]

    original_element_i = data["element_i"].copy()
    original_element_j = data["element_j"].copy()

    data.loc[swap_mask, "element_i"] = original_element_j[swap_mask]
    data.loc[swap_mask, "element_j"] = original_element_i[swap_mask]

    data = data.sort_values(
        ["system", "pair_i", "pair_j"]
    ).drop_duplicates(
        subset=["system", "pair_i", "pair_j"],
        keep="first",
    )

    return data


def load_mapping(path: Path) -> tuple[dict[int, int], dict[int, int]]:
    """Load parent-to-functionalized and reverse atom mappings."""

    mapping = pd.read_csv(path)

    required = {"parent_atom", "functionalized_atom"}
    missing = required - set(mapping.columns)

    if missing:
        raise ValueError(
            "Atom mapping is missing required columns: "
            + ", ".join(sorted(missing))
        )

    matched = mapping.dropna(
        subset=["parent_atom", "functionalized_atom"]
    ).copy()

    matched["parent_atom"] = matched["parent_atom"].astype(int)
    matched["functionalized_atom"] = matched[
        "functionalized_atom"
    ].astype(int)

    parent_to_functionalized = dict(
        zip(
            matched["parent_atom"],
            matched["functionalized_atom"],
            strict=True,
        )
    )

    functionalized_to_parent = {
        functionalized: parent
        for parent, functionalized
        in parent_to_functionalized.items()
    }

    return parent_to_functionalized, functionalized_to_parent


def make_lookup(
    system_inventory: pd.DataFrame,
) -> dict[tuple[int, int], dict[str, object]]:
    """Index a system's bond records by canonical atom pair."""

    return {
        (int(row["pair_i"]), int(row["pair_j"])): row
        for row in system_inventory.to_dict("records")
    }


def generate_rows(
    inventory: pd.DataFrame,
    parent_to_functionalized: dict[int, int],
    functionalized_to_parent: dict[int, int],
    parent_system: str,
    functionalized_system: str,
) -> list[dict[str, object]]:
    """Generate matched and functionalized-only candidate bonds."""

    parent = inventory[
        inventory["system"] == parent_system
    ].copy()

    functionalized = inventory[
        inventory["system"] == functionalized_system
    ].copy()

    if parent.empty:
        raise ValueError(
            f"No bond records found for parent system {parent_system!r}."
        )

    if functionalized.empty:
        raise ValueError(
            "No bond records found for functionalized system "
            f"{functionalized_system!r}."
        )

    parent_lookup = make_lookup(parent)
    functionalized_lookup = make_lookup(functionalized)

    rows: list[dict[str, object]] = []
    counters: dict[str, int] = {}

    def next_id(prefix: str) -> str:
        counters[prefix] = counters.get(prefix, 0) + 1
        return f"{prefix.lower()}_{counters[prefix]}"

    # Bonds that can be compared directly between the two systems.
    for (parent_i, parent_j), parent_row in sorted(
        parent_lookup.items()
    ):
        if (
            parent_i not in parent_to_functionalized
            or parent_j not in parent_to_functionalized
        ):
            continue

        functionalized_pair = canonical_pair(
            parent_to_functionalized[parent_i],
            parent_to_functionalized[parent_j],
        )

        if functionalized_pair not in functionalized_lookup:
            continue

        bond_label = element_pair(
            parent_row["element_i"],
            parent_row["element_j"],
        )

        rows.append(
            {
                "bond_id": next_id(bond_label.replace("-", "")),
                "label": bond_label,
                "parent_i": parent_i,
                "parent_j": parent_j,
                "functionalized_i": functionalized_pair[0],
                "functionalized_j": functionalized_pair[1],
                "comparison_type": "matched",
            }
        )

    mapped_functionalized_atoms = set(
        functionalized_to_parent
    )

    # Bonds involving at least one atom unique to the functionalized molecule.
    for (
        functionalized_i,
        functionalized_j,
    ), functionalized_row in sorted(functionalized_lookup.items()):
        both_atoms_mapped = (
            functionalized_i in mapped_functionalized_atoms
            and functionalized_j in mapped_functionalized_atoms
        )

        if both_atoms_mapped:
            continue

        bond_label = element_pair(
            functionalized_row["element_i"],
            functionalized_row["element_j"],
        )

        rows.append(
            {
                "bond_id": next_id(
                    f"new{bond_label.replace('-', '')}"
                ),
                "label": f"new {bond_label}",
                "parent_i": pd.NA,
                "parent_j": pd.NA,
                "functionalized_i": functionalized_i,
                "functionalized_j": functionalized_j,
                "comparison_type": "functionalized_only",
            }
        )

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a candidate selected_bonds.csv from the parsed "
            "Mayer bond inventory and atom mapping."
        )
    )

    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path(
            "results/ground_state/bond_orders/bond_inventory.csv"
        ),
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("config/pdi_core_atom_mapping.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("config/selected_bonds_candidates.csv"),
    )
    parser.add_argument(
        "--parent-system",
        default="pdi",
    )
    parser.add_argument(
        "--functionalized-system",
        default="pdi_terminal_functionalized",
    )
    parser.add_argument(
        "--minimum-bond-order",
        type=float,
        default=0.50,
        help=(
            "Minimum positive Mayer bond order used to identify "
            "candidate covalent bonds."
        ),
    )
    parser.add_argument(
        "--include-hydrogen",
        action="store_true",
        help="Include bonds involving hydrogen atoms.",
    )

    args = parser.parse_args()

    if not args.inventory.exists():
        raise FileNotFoundError(
            f"Bond inventory not found: {args.inventory}"
        )

    if not args.mapping.exists():
        raise FileNotFoundError(
            f"Atom mapping not found: {args.mapping}"
        )

    inventory = pd.read_csv(args.inventory)

    inventory = prepare_inventory(
        inventory=inventory,
        minimum_bond_order=args.minimum_bond_order,
        include_hydrogen=args.include_hydrogen,
    )

    (
        parent_to_functionalized,
        functionalized_to_parent,
    ) = load_mapping(args.mapping)

    rows = generate_rows(
        inventory=inventory,
        parent_to_functionalized=parent_to_functionalized,
        functionalized_to_parent=functionalized_to_parent,
        parent_system=args.parent_system,
        functionalized_system=args.functionalized_system,
    )

    output = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    matched_count = (
        output["comparison_type"].eq("matched").sum()
        if not output.empty
        else 0
    )

    functionalized_only_count = (
        output["comparison_type"]
        .eq("functionalized_only")
        .sum()
        if not output.empty
        else 0
    )

    print(f"Wrote {len(output)} candidate bonds to {args.output}")
    print(f"Matched bonds: {matched_count}")
    print(
        "Functionalized-only bonds: "
        f"{functionalized_only_count}"
    )
    print(
        "Review this candidate file manually before replacing "
        "config/selected_bonds.csv."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
