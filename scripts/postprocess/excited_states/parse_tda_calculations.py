from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

try:
    from scripts.postprocess.ground_state.common import parse_float, read_csv_dicts, read_text, repo_root, write_csv
except ImportError:
    from common import parse_float, read_csv_dicts, read_text, repo_root, write_csv

EV_NM = 1239.841984

SYSTEM_LABELS = {
    "pdi": "Parent PDI",
    "pdi_terminal_functionalized": "Terminal-functionalized PDI",
}

MULTIPLICITY_LABELS = {
    "singlet": "Singlet",
    "triplet": "Triplet",
}

MULTIPLICITY_NUMBERS = {
    "singlet": 1,
    "triplet": 3,
}

SECTION_NAMES = {
    "singlet": "SINGLETS",
    "triplet": "TRIPLETS",
}

STATE_PREFIXES = {
    "singlet": "S",
    "triplet": "T",
}

TDA_OUTPUTS = {
    "pdi": {
        "singlet": Path("calculations/pdi/excited_state_calculations/tda_singlets/pdi_tda_singlets_tprint.out"),
        "triplet": Path("calculations/pdi/excited_state_calculations/tda_triplets/pdi_tda_triplets_tprint.out"),
    },
    "pdi_terminal_functionalized": {
        "singlet": Path(
            "calculations/pdi_terminal_functionalized/excited_state_calculations/tda_singlets/"
            "pdi_terminal_functionalized_tda_singlets_tprint.out"
        ),
        "triplet": Path(
            "calculations/pdi_terminal_functionalized/excited_state_calculations/tda_triplets/"
            "pdi_terminal_functionalized_tda_triplets_tprint.out"
        ),
    },
}

STATE_RE = re.compile(
    r"^STATE\s+(\d+):\s+E=\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+au\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+eV\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+cm\*\*-1\s+"
    r"<S\*\*2>\s+=\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"Mult\s+(\d+)",
    re.MULTILINE,
)

TRANSITION_RE = re.compile(
    r"^\s*(\d+)([abAB]?)\s*->\s*(\d+)([abAB]?)\s*:\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"\(c=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\)",
    re.MULTILINE,
)

ABSORPTION_RE = re.compile(
    r"^\s*\S+\s+->\s+(\d+)-\S+\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)",
    re.MULTILINE,
)


class TDAParseError(ValueError):
    """Raised when a TDA output cannot be parsed safely."""


def relative_path(path: Path, repo: Path) -> str:
    return str(path.relative_to(repo))


def validate_normal_termination(text: str, system: str, multiplicity: str, path: Path) -> None:
    if "****ORCA TERMINATED NORMALLY****" not in text:
        raise TDAParseError(f"{system} {multiplicity} output did not terminate normally: {path}")


def extract_state_section(text: str, system: str, multiplicity: str, path: Path) -> tuple[str, int | None]:
    target = SECTION_NAMES[multiplicity]
    header_re = re.compile(r"TD-DFT/TDA EXCITED STATES \((SINGLETS|TRIPLETS)\)")
    headers = list(header_re.finditer(text))
    selected = [match for match in headers if match.group(1) == target]
    if not selected:
        raise TDAParseError(f"No {target} excited-state section for {system} {multiplicity}: {path}")
    if len(selected) > 1:
        raise TDAParseError(f"Multiple {target} excited-state sections for {system} {multiplicity}: {path}")

    match = selected[0]
    start = match.end()
    following_headers = [other.start() for other in headers if other.start() > start]
    spectra = text.find("TD-DFT/TDA-EXCITATION SPECTRA", start)
    ends = following_headers + ([spectra] if spectra != -1 else [])
    end = min(ends) if ends else len(text)

    roots = re.findall(r"Number of roots to be determined\s+\.\.\.\s+(\d+)", text[: match.start()])
    requested_roots = int(roots[-1]) if roots else None
    return text[start:end], requested_roots


def parse_absorption_section(text: str) -> dict[int, dict[str, float]]:
    marker = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    start = text.find(marker)
    if start == -1:
        return {}
    next_marker = text.find("ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS", start)
    section = text[start: next_marker if next_marker != -1 else len(text)]
    rows = {}
    for match in ABSORPTION_RE.finditer(section):
        state_index = int(match.group(1))
        rows[state_index] = {
            "excitation_energy_ev": parse_float(match.group(2)),
            "excitation_energy_cm1": parse_float(match.group(3)),
            "wavelength_nm": parse_float(match.group(4)),
            "oscillator_strength": parse_float(match.group(5)),
            "transition_dipole_strength_au2": parse_float(match.group(6)),
            "transition_dipole_x_au": parse_float(match.group(7)),
            "transition_dipole_y_au": parse_float(match.group(8)),
            "transition_dipole_z_au": parse_float(match.group(9)),
        }
    return rows


def read_molden_homo_indices(repo: Path) -> dict[str, int]:
    path = repo / "results/ground_state/orbitals/frontier_orbital_energies.csv"
    if not path.is_file():
        return {}
    indices = {}
    for row in read_csv_dicts(path):
        system = row.get("system", "").strip()
        homo = row.get("homo_index", "").strip()
        if system and homo:
            indices[system] = int(float(homo))
    return indices


def transition_homo_index(
    text: str,
    system: str,
    multiplicity: str,
    path: Path,
    molden_homo_indices: dict[str, int],
) -> int:
    match = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", text)
    if match:
        electrons = int(match.group(1))
        if electrons % 2:
            raise TDAParseError(f"Odd electron count cannot define a restricted HOMO index in {path}: {electrons}")
        return electrons // 2 - 1

    if system in molden_homo_indices:
        # ORCA TDA transition listings in these outputs are zero-based, while
        # Molden orbital tables are one-based.
        return molden_homo_indices[system] - 1

    raise TDAParseError(
        f"Could not determine ORCA transition HOMO index for {system} {multiplicity}; "
        f"missing NEL in {path} and frontier_orbital_energies.csv entry"
    )


def orbital_label(orbital_index: int, homo_index: int) -> str:
    if orbital_index <= homo_index:
        offset = homo_index - orbital_index
        return "HOMO" if offset == 0 else f"HOMO-{offset}"
    offset = orbital_index - homo_index - 1
    return "LUMO" if offset == 0 else f"LUMO+{offset}"


def parse_transition_lines(
    block: str,
    system: str,
    multiplicity: str,
    state_index: int,
    state_label: str,
    source_file: str,
    homo_index: int,
) -> list[dict[str, Any]]:
    rows = []
    for match in TRANSITION_RE.finditer(block):
        occupied = int(match.group(1))
        virtual = int(match.group(3))
        signed_coefficient = parse_float(match.group(6))
        weight = parse_float(match.group(5))
        occupied_label = orbital_label(occupied, homo_index)
        virtual_label = orbital_label(virtual, homo_index)
        rows.append(
            {
                "system": system,
                "multiplicity": multiplicity,
                "state_index": state_index,
                "state_label": state_label,
                "source_file": source_file,
                "occupied_orbital_index": occupied,
                "virtual_orbital_index": virtual,
                "occupied_spin": match.group(2).lower(),
                "virtual_spin": match.group(4).lower(),
                "occupied_orbital_label": occupied_label,
                "virtual_orbital_label": virtual_label,
                "transition_label": f"{occupied_label} → {virtual_label}",
                "transition_coefficient": abs(signed_coefficient),
                "signed_coefficient": signed_coefficient,
                "squared_coefficient": signed_coefficient * signed_coefficient,
                "transition_weight": weight,
            }
        )
    return rows


def dominant_transition_summary(transitions: list[dict[str, Any]], minimum_weight: float) -> str:
    kept = [row for row in transitions if row["transition_weight"] >= minimum_weight]
    kept.sort(key=lambda row: row["transition_weight"], reverse=True)
    return "; ".join(f"{row['transition_label']} ({100 * row['transition_weight']:.1f}%)" for row in kept)


def parse_state_section(
    section: str,
    system: str,
    multiplicity: str,
    source_file: str,
    homo_index: int,
    minimum_transition_weight: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    matches = list(STATE_RE.finditer(section))
    if not matches:
        raise TDAParseError(f"No STATE blocks parsed for {system} {multiplicity}: {source_file}")

    states = []
    transitions = []
    expected_mult = MULTIPLICITY_NUMBERS[multiplicity]
    prefix = STATE_PREFIXES[multiplicity]

    for i, match in enumerate(matches, 1):
        state_block_end = matches[i].start() if i < len(matches) else len(section)
        transition_block = section[match.end():state_block_end]
        orca_state_index = int(match.group(1))
        parsed_mult = int(match.group(6))
        if parsed_mult != expected_mult:
            raise TDAParseError(
                f"Unexpected multiplicity in {system} {multiplicity} state block: "
                f"STATE {orca_state_index} has Mult {parsed_mult}"
            )

        state_label = f"{prefix}{i}"
        state_transitions = parse_transition_lines(
            transition_block,
            system,
            multiplicity,
            i,
            state_label,
            source_file,
            homo_index,
        )
        transitions.extend(state_transitions)
        energy_ev = parse_float(match.group(3))
        states.append(
            {
                "system": system,
                "system_label": SYSTEM_LABELS[system],
                "multiplicity": multiplicity,
                "multiplicity_label": MULTIPLICITY_LABELS[multiplicity],
                "state_index": i,
                "state_label": state_label,
                "orca_state_index": orca_state_index,
                "source_file": source_file,
                "excitation_energy_au": parse_float(match.group(2)),
                "excitation_energy_ev": energy_ev,
                "excitation_energy_cm1": parse_float(match.group(4)),
                "wavelength_nm": EV_NM / energy_ev,
                "spin_expectation_s2": parse_float(match.group(5)),
                "oscillator_strength": math.nan,
                "transition_dipole_strength_au2": math.nan,
                "transition_dipole_x_au": math.nan,
                "transition_dipole_y_au": math.nan,
                "transition_dipole_z_au": math.nan,
                "dominant_transitions": dominant_transition_summary(state_transitions, minimum_transition_weight),
            }
        )
    return states, transitions


def apply_oscillator_strengths(
    states: list[dict[str, Any]],
    absorption_rows: dict[int, dict[str, float]],
    system: str,
    multiplicity: str,
    path: Path,
) -> list[str]:
    warnings = []
    if multiplicity == "triplet":
        # Triplet excitations are spin-forbidden in this TDA workflow, so keep
        # oscillator strengths as NaN rather than inventing zero-intensity data.
        return warnings

    if not absorption_rows:
        raise TDAParseError(f"No electric-dipole absorption spectrum parsed for {system} singlet: {path}")

    if len(absorption_rows) != len(states):
        raise TDAParseError(
            f"Singlet oscillator-strength row count mismatch for {system}: "
            f"{len(absorption_rows)} spectrum rows but {len(states)} states in {path}"
        )

    for state in states:
        state_index = state["state_index"]
        spectrum = absorption_rows.get(state_index)
        if spectrum is None:
            raise TDAParseError(f"Missing oscillator strength for {system} S{state_index}: {path}")
        state.update(spectrum)
        state["wavelength_nm"] = EV_NM / state["excitation_energy_ev"]
    return warnings


def validate_rows(
    states: list[dict[str, Any]],
    transitions: list[dict[str, Any]],
    system: str,
    multiplicity: str,
    path: Path,
) -> None:
    state_indices = [row["state_index"] for row in states]
    if len(state_indices) != len(set(state_indices)):
        raise TDAParseError(f"Duplicate state indices for {system} {multiplicity}: {path}")

    state_keys = {(row["system"], row["multiplicity"], row["state_index"]) for row in states}
    for row in states:
        energy = row["excitation_energy_ev"]
        wavelength = row["wavelength_nm"]
        if not math.isfinite(energy) or energy <= 0:
            raise TDAParseError(f"Invalid excitation energy for {row['state_label']} in {path}: {energy}")
        if abs(wavelength - EV_NM / energy) > 1e-8:
            raise TDAParseError(f"Wavelength mismatch for {row['state_label']} in {path}")
        fosc = row["oscillator_strength"]
        if multiplicity == "singlet" and (not math.isfinite(fosc) or fosc < 0):
            raise TDAParseError(f"Invalid singlet oscillator strength for {row['state_label']} in {path}: {fosc}")

    for row in transitions:
        key = (row["system"], row["multiplicity"], row["state_index"])
        if key not in state_keys:
            raise TDAParseError(f"Transition row maps to missing state {key}: {path}")
        if row["transition_weight"] < 0 or not math.isfinite(row["transition_weight"]):
            raise TDAParseError(f"Invalid transition weight for {key}: {row['transition_weight']}")


def parse_output(
    repo: Path,
    system: str,
    multiplicity: str,
    path: Path,
    molden_homo_indices: dict[str, int],
    minimum_transition_weight: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    absolute_path = repo / path
    text = read_text(absolute_path)
    validate_normal_termination(text, system, multiplicity, absolute_path)
    section, requested_roots = extract_state_section(text, system, multiplicity, absolute_path)
    source_file = relative_path(absolute_path, repo)
    homo_index = transition_homo_index(text, system, multiplicity, absolute_path, molden_homo_indices)

    states, transitions = parse_state_section(
        section,
        system,
        multiplicity,
        source_file,
        homo_index,
        minimum_transition_weight,
    )
    warnings = apply_oscillator_strengths(
        states,
        parse_absorption_section(text),
        system,
        multiplicity,
        absolute_path,
    )
    validate_rows(states, transitions, system, multiplicity, absolute_path)

    finite_oscillator = [
        row["oscillator_strength"] for row in states if isinstance(row["oscillator_strength"], float) and math.isfinite(row["oscillator_strength"])
    ]
    summary = {
        "system": system,
        "system_label": SYSTEM_LABELS[system],
        "multiplicity": multiplicity,
        "multiplicity_label": MULTIPLICITY_LABELS[multiplicity],
        "source_output_path": source_file,
        "normal_termination": True,
        "requested_or_detected_number_of_roots": requested_roots if requested_roots is not None else "",
        "parsed_number_of_states": len(states),
        "minimum_excitation_energy_ev": min(row["excitation_energy_ev"] for row in states),
        "maximum_excitation_energy_ev": max(row["excitation_energy_ev"] for row in states),
        "number_of_states_with_oscillator_strengths": len(finite_oscillator),
        "parser_warnings": "; ".join(warnings),
    }
    return states, transitions, summary


def export_rows(repo: Path, states: list[dict[str, Any]], transitions: list[dict[str, Any]], summary: list[dict[str, Any]]) -> None:
    out_dir = repo / "results/excited_state/tda_uv_vis"
    write_csv(
        out_dir / "excitation_states.csv",
        states,
        [
            "system",
            "system_label",
            "multiplicity",
            "multiplicity_label",
            "state_index",
            "state_label",
            "orca_state_index",
            "source_file",
            "excitation_energy_au",
            "excitation_energy_ev",
            "excitation_energy_cm1",
            "wavelength_nm",
            "spin_expectation_s2",
            "oscillator_strength",
            "transition_dipole_strength_au2",
            "transition_dipole_x_au",
            "transition_dipole_y_au",
            "transition_dipole_z_au",
            "dominant_transitions",
        ],
    )
    write_csv(
        out_dir / "excitation_transitions.csv",
        transitions,
        [
            "system",
            "multiplicity",
            "state_index",
            "state_label",
            "source_file",
            "occupied_orbital_index",
            "virtual_orbital_index",
            "occupied_spin",
            "virtual_spin",
            "occupied_orbital_label",
            "virtual_orbital_label",
            "transition_label",
            "transition_coefficient",
            "signed_coefficient",
            "squared_coefficient",
            "transition_weight",
        ],
    )
    write_csv(
        out_dir / "tda_parse_summary.csv",
        summary,
        [
            "system",
            "system_label",
            "multiplicity",
            "multiplicity_label",
            "source_output_path",
            "normal_termination",
            "requested_or_detected_number_of_roots",
            "parsed_number_of_states",
            "minimum_excitation_energy_ev",
            "maximum_excitation_energy_ev",
            "number_of_states_with_oscillator_strengths",
            "parser_warnings",
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--minimum-transition-weight", type=float, default=0.05)
    args = parser.parse_args()

    repo = repo_root(args.repo)
    molden_homo_indices = read_molden_homo_indices(repo)
    all_states: list[dict[str, Any]] = []
    all_transitions: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    for system, paths in TDA_OUTPUTS.items():
        for multiplicity, path in paths.items():
            states, transitions, summary = parse_output(
                repo,
                system,
                multiplicity,
                path,
                molden_homo_indices,
                args.minimum_transition_weight,
            )
            all_states.extend(states)
            all_transitions.extend(transitions)
            summaries.append(summary)
            print(
                f"{system} {multiplicity}: parsed {len(states)} states, "
                f"{len(transitions)} transitions from {path}"
            )

    export_rows(repo, all_states, all_transitions, summaries)
    print(f"Wrote {repo / 'results/excited_state/tda_uv_vis/excitation_states.csv'}")
    print(f"Wrote {repo / 'results/excited_state/tda_uv_vis/excitation_transitions.csv'}")
    print(f"Wrote {repo / 'results/excited_state/tda_uv_vis/tda_parse_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
