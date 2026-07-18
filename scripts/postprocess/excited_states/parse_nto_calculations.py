from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

try:
    from scripts.postprocess.ground_state.common import parse_float, repo_root, write_csv, write_json
except ImportError:
    from common import parse_float, repo_root, write_csv, write_json

SYSTEM_LABELS = {
    "pdi": "Parent PDI",
    "pdi_terminal_functionalized": "Terminal-functionalized PDI",
}

MULTIPLICITY_LABELS = {
    "singlet": "Singlet",
    "triplet": "Triplet",
}

STATE_PREFIXES = {
    "singlet": "S",
    "triplet": "T",
}

DEDICATED_MULTIPLICITY_PATTERNS = {
    "singlet": ("nto_singlet", "singlet"),
    "triplet": ("nto_triplet", "triplet"),
}

FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"

STATE_RE = re.compile(
    rf"^STATE\s+(\d+):\s+E=\s+({FLOAT_RE})\s+au\s+({FLOAT_RE})\s+eV\s+"
    rf"({FLOAT_RE})\s+cm\*\*-1(?:\s+<S\*\*2>\s+=\s+({FLOAT_RE}))?(?:\s+Mult\s+(\d+))?",
    re.MULTILINE,
)

SECTION_RE = re.compile(r"TD-DFT/TDA EXCITED STATES \((SINGLETS|TRIPLETS)\)")
NTO_BLOCK_RE = re.compile(r"NATURAL TRANSITION ORBITALS FOR STATE\s+(\d+)")
NTO_FILE_RE = re.compile(r"Natural Transition Orbitals were saved in\s+(\S+)")
NTO_THRESHOLD_RE = re.compile(r"Threshold for printing occupation numbers\s+(" + FLOAT_RE + r")")
NTO_ENERGY_RE = re.compile(
    rf"^\s*E=\s+({FLOAT_RE})\s+au\s+({FLOAT_RE})\s+eV\s+({FLOAT_RE})\s+cm\*\*-1",
    re.MULTILINE,
)
NTO_PAIR_RE = re.compile(rf"^\s*(\d+[a-zA-Z]?)\s*->\s*(\d+[a-zA-Z]?)\s+:\s+n=\s+({FLOAT_RE})", re.MULTILINE)


class NTOParseError(ValueError):
    """Raised when an NTO output contains a critical validation error."""


@dataclass(frozen=True)
class StateRecord:
    system: str
    source_output: Path
    global_state_number: int
    multiplicity: str
    local_state_index: int
    state_label: str
    energy_hartree: float
    energy_ev: float
    energy_cm1: float
    spin_expectation_s2: float | None
    multiplicity_value: int | None


@dataclass
class NTOPair:
    pair_index: int
    hole_orbital: str
    electron_orbital: str
    weight: float
    cumulative_weight: float


@dataclass
class NTOOccurrence:
    system: str
    source_output: Path
    dedicated_multiplicity: str
    job_normal_termination: bool
    global_state_number: int
    state: StateRecord | None
    nto_block_energy_hartree: float | None
    nto_block_energy_ev: float | None
    nto_block_energy_cm1: float | None
    nto_filename: str
    nto_path: Path
    nto_file_exists: bool
    nto_print_threshold: float | None
    pairs: list[NTOPair]
    warnings: list[str] = field(default_factory=list)
    duplicate_status: str = "unique"
    discarded_duplicate_sources: list[str] = field(default_factory=list)


def relative_path(path: Path, repo: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo.resolve()))
    except ValueError:
        return str(path)


def resolve_cli_path(repo: Path, value: Path | None, default: Path) -> Path:
    path = value or default
    return path if path.is_absolute() else repo / path


def infer_system(path: Path) -> str | None:
    normalized = path.as_posix().lower().replace("__", "_")
    if "pdi_terminal_functionalized" in normalized:
        return "pdi_terminal_functionalized"
    if re.search(r"(^|/)pdi(/|_|-)", normalized) or "pdi_nto" in normalized:
        return "pdi"
    return None


def infer_dedicated_multiplicity(path: Path) -> str:
    normalized = path.as_posix().lower()
    for multiplicity, patterns in DEDICATED_MULTIPLICITY_PATTERNS.items():
        if any(pattern in normalized for pattern in patterns):
            return multiplicity
    return ""


def discover_nto_outputs(calculation_root: Path) -> list[Path]:
    outputs = []
    for path in sorted(calculation_root.rglob("*.out")):
        normalized = path.as_posix().lower()
        if "nto" not in normalized:
            continue
        if infer_system(path) is None:
            continue
        outputs.append(path)
    return outputs


def parse_state_sections(text: str, system: str, source_output: Path) -> dict[int, StateRecord]:
    states: dict[int, StateRecord] = {}
    headers = list(SECTION_RE.finditer(text))
    nto_start = text.find("NATURAL TRANSITION ORBITALS")

    for header_index, header in enumerate(headers):
        section_name = header.group(1)
        multiplicity = "singlet" if section_name == "SINGLETS" else "triplet"
        prefix = STATE_PREFIXES[multiplicity]
        next_header_start = headers[header_index + 1].start() if header_index + 1 < len(headers) else len(text)
        end_candidates = [next_header_start]
        if nto_start != -1 and nto_start > header.end():
            end_candidates.append(nto_start)
        section = text[header.end() : min(end_candidates)]

        for local_index, match in enumerate(STATE_RE.finditer(section), 1):
            global_state_number = int(match.group(1))
            states[global_state_number] = StateRecord(
                system=system,
                source_output=source_output,
                global_state_number=global_state_number,
                multiplicity=multiplicity,
                local_state_index=local_index,
                state_label=f"{prefix}{local_index}",
                energy_hartree=parse_float(match.group(2)),
                energy_ev=parse_float(match.group(3)),
                energy_cm1=parse_float(match.group(4)),
                spin_expectation_s2=parse_float(match.group(5)) if match.group(5) else None,
                multiplicity_value=int(match.group(6)) if match.group(6) else None,
            )
    return states


def parse_nto_blocks(
    text: str,
    repo: Path,
    output_path: Path,
    system: str,
    states: dict[int, StateRecord],
    energy_tolerance_ev: float,
) -> list[NTOOccurrence]:
    blocks = list(NTO_BLOCK_RE.finditer(text))
    occurrences = []
    dedicated_multiplicity = infer_dedicated_multiplicity(output_path)
    normal_termination = "****ORCA TERMINATED NORMALLY****" in text

    for block_index, header in enumerate(blocks):
        block = text[header.start() : blocks[block_index + 1].start() if block_index + 1 < len(blocks) else len(text)]
        global_state_number = int(header.group(1))
        state = states.get(global_state_number)
        warnings: list[str] = []

        filename_match = NTO_FILE_RE.search(block)
        if filename_match is None:
            nto_filename = ""
            nto_path = output_path.parent
            warnings.append("missing_nto_filename")
        else:
            nto_filename = filename_match.group(1)
            nto_path = output_path.parent / nto_filename

        threshold_match = NTO_THRESHOLD_RE.search(block)
        energy_match = NTO_ENERGY_RE.search(block)
        nto_block_energy_hartree = parse_float(energy_match.group(1)) if energy_match else None
        nto_block_energy_ev = parse_float(energy_match.group(2)) if energy_match else None
        nto_block_energy_cm1 = parse_float(energy_match.group(3)) if energy_match else None

        if state is None:
            warnings.append("unknown_state_number")
        elif nto_block_energy_ev is not None and abs(nto_block_energy_ev - state.energy_ev) > energy_tolerance_ev:
            warnings.append("energy_mismatch")

        pairs = []
        cumulative_weight = 0.0
        for pair_index, pair_match in enumerate(NTO_PAIR_RE.finditer(block), 1):
            weight = parse_float(pair_match.group(3))
            cumulative_weight += weight
            pairs.append(
                NTOPair(
                    pair_index=pair_index,
                    hole_orbital=pair_match.group(1),
                    electron_orbital=pair_match.group(2),
                    weight=weight,
                    cumulative_weight=cumulative_weight,
                )
            )

        if not pairs:
            warnings.append("no_nto_pairs")
        if not normal_termination:
            warnings.append("incomplete_orca_job")
        if nto_filename and not nto_path.is_file():
            warnings.append("missing_nto_file")
        if state is not None and dedicated_multiplicity and dedicated_multiplicity != state.multiplicity:
            warnings.append("source_dedicated_multiplicity_mismatch")

        occurrences.append(
            NTOOccurrence(
                system=system,
                source_output=output_path,
                dedicated_multiplicity=dedicated_multiplicity,
                job_normal_termination=normal_termination,
                global_state_number=global_state_number,
                state=state,
                nto_block_energy_hartree=nto_block_energy_hartree,
                nto_block_energy_ev=nto_block_energy_ev,
                nto_block_energy_cm1=nto_block_energy_cm1,
                nto_filename=nto_filename,
                nto_path=nto_path,
                nto_file_exists=nto_path.is_file() if nto_filename else False,
                nto_print_threshold=parse_float(threshold_match.group(1)) if threshold_match else None,
                pairs=pairs,
                warnings=warnings,
            )
        )

    return occurrences


def pair_signature_conflicts(left: NTOOccurrence, right: NTOOccurrence, weight_tolerance: float = 1e-4) -> bool:
    if not left.pairs or not right.pairs:
        return False
    if len(left.pairs) != len(right.pairs):
        return True
    for left_pair, right_pair in zip(left.pairs, right.pairs, strict=True):
        if left_pair.hole_orbital != right_pair.hole_orbital:
            return True
        if left_pair.electron_orbital != right_pair.electron_orbital:
            return True
        if abs(left_pair.weight - right_pair.weight) > weight_tolerance:
            return True
    return False


def occurrence_preference_key(repo: Path, occurrence: NTOOccurrence) -> tuple[int, int, int, str]:
    dedicated_match = int(bool(occurrence.state and occurrence.dedicated_multiplicity == occurrence.state.multiplicity))
    return (
        int(occurrence.job_normal_termination),
        int(occurrence.nto_file_exists),
        dedicated_match,
        relative_path(occurrence.source_output, repo),
    )


def canonical_identity(occurrence: NTOOccurrence) -> tuple[str, str, int] | None:
    if occurrence.state is None:
        return None
    return (
        occurrence.system,
        occurrence.state.multiplicity,
        occurrence.state.local_state_index,
    )


def deduplicate_occurrences(
    repo: Path,
    occurrences: Iterable[NTOOccurrence],
    energy_tolerance_ev: float,
) -> tuple[list[NTOOccurrence], dict[str, str]]:
    groups: dict[tuple[str, str, int], list[NTOOccurrence]] = {}
    for occurrence in occurrences:
        identity = canonical_identity(occurrence)
        if identity is None:
            continue
        groups.setdefault(identity, []).append(occurrence)

    canonical = []
    source_status_values: dict[str, set[str]] = {}

    def add_source_status(path: Path, status: str) -> None:
        source_status_values.setdefault(relative_path(path, repo), set()).add(status)

    for group in groups.values():
        if len(group) == 1:
            selected = group[0]
            selected.duplicate_status = "unique"
            canonical.append(selected)
            add_source_status(selected.source_output, "unique")
            continue

        selected = sorted(group, key=lambda item: occurrence_preference_key(repo, item), reverse=True)[0]
        energies = [item.state.energy_ev for item in group if item.state is not None]
        conflict = max(energies) - min(energies) > energy_tolerance_ev
        conflict = conflict or any(pair_signature_conflicts(selected, item) for item in group if item is not selected)

        duplicate_sources = [
            relative_path(item.source_output, repo)
            for item in group
            if item is not selected
        ]
        selected.discarded_duplicate_sources = duplicate_sources
        selected.duplicate_status = "duplicate_conflict_selected" if conflict else "duplicate_selected"
        if conflict and "duplicate_conflict" not in selected.warnings:
            selected.warnings.append("duplicate_conflict")
        elif "duplicate_state" not in selected.warnings:
            selected.warnings.append("duplicate_state")
        canonical.append(selected)

        add_source_status(selected.source_output, selected.duplicate_status)
        discarded_status = "duplicate_conflict_discarded" if conflict else "duplicate_discarded"
        for item in group:
            if item is selected:
                continue
            item.duplicate_status = discarded_status
            add_source_status(item.source_output, discarded_status)

    return sorted(
        canonical,
        key=lambda item: (
            item.system,
            item.state.multiplicity if item.state else "",
            item.state.local_state_index if item.state else 0,
        ),
    ), {source: "; ".join(sorted(statuses)) for source, statuses in source_status_values.items()}


def minimum_pairs_for_cutoff(pairs: list[NTOPair], cutoff: float) -> tuple[int, float, str]:
    if not pairs:
        return 0, math.nan, "no_pairs"
    for pair in pairs:
        if pair.cumulative_weight >= cutoff:
            return pair.pair_index, pair.cumulative_weight, "reached"
    return len(pairs), pairs[-1].cumulative_weight, "printed_weights_below_cutoff"


def pair_indices(count: int) -> str:
    return ";".join(str(index) for index in range(1, count + 1))


def classify_nto_character(leading_pair_weight: float) -> str:
    if leading_pair_weight >= 0.90:
        return "essentially single-pair"
    if leading_pair_weight >= 0.75:
        return "dominant primary pair with secondary mixing"
    return "strongly mixed multi-pair state"


def warning_text(warnings: Iterable[str]) -> str:
    return "; ".join(sorted(set(warning for warning in warnings if warning)))


def occurrence_warning_values(occurrence: NTOOccurrence) -> list[str]:
    values = list(occurrence.warnings)
    if occurrence.duplicate_status != "unique":
        values.append(occurrence.duplicate_status)
    return values


def raw_pair_rows(repo: Path, occurrences: Iterable[NTOOccurrence]) -> list[dict[str, Any]]:
    rows = []
    for occurrence in occurrences:
        state = occurrence.state
        for pair in occurrence.pairs:
            rows.append(
                {
                    "system": occurrence.system,
                    "system_label": SYSTEM_LABELS[occurrence.system],
                    "source_output": relative_path(occurrence.source_output, repo),
                    "job_normal_termination": occurrence.job_normal_termination,
                    "multiplicity": state.multiplicity if state else "",
                    "multiplicity_label": MULTIPLICITY_LABELS[state.multiplicity] if state else "",
                    "global_state_number": occurrence.global_state_number,
                    "local_state_index": state.local_state_index if state else "",
                    "state_label": state.state_label if state else "",
                    "state_energy_hartree": state.energy_hartree if state else "",
                    "state_energy_ev": state.energy_ev if state else "",
                    "state_energy_cm1": state.energy_cm1 if state else "",
                    "nto_block_energy_ev": occurrence.nto_block_energy_ev if occurrence.nto_block_energy_ev is not None else "",
                    "nto_filename": occurrence.nto_filename,
                    "nto_path": relative_path(occurrence.nto_path, repo) if occurrence.nto_filename else "",
                    "nto_file_exists": occurrence.nto_file_exists,
                    "nto_print_threshold": occurrence.nto_print_threshold if occurrence.nto_print_threshold is not None else "",
                    "pair_index": pair.pair_index,
                    "hole_orbital": pair.hole_orbital,
                    "electron_orbital": pair.electron_orbital,
                    "weight": pair.weight,
                    "cumulative_weight": pair.cumulative_weight,
                    "parse_warning": warning_text(occurrence_warning_values(occurrence)),
                }
            )
    return rows


def state_summary_rows(
    repo: Path,
    canonical: Iterable[NTOOccurrence],
    main_cutoff: float,
    si_cutoff: float,
) -> list[dict[str, Any]]:
    rows = []
    for occurrence in canonical:
        if occurrence.state is None:
            continue
        pairs = occurrence.pairs
        leading_weight = pairs[0].weight if pairs else math.nan
        total_weight = pairs[-1].cumulative_weight if pairs else math.nan
        main_count, main_weight, main_status = minimum_pairs_for_cutoff(pairs, main_cutoff)
        si_count, si_weight, si_status = minimum_pairs_for_cutoff(pairs, si_cutoff)
        status = "ok" if main_status == "reached" and si_status == "reached" else warning_text([main_status, si_status])
        rows.append(
            {
                "system": occurrence.system,
                "system_label": SYSTEM_LABELS[occurrence.system],
                "multiplicity": occurrence.state.multiplicity,
                "multiplicity_label": MULTIPLICITY_LABELS[occurrence.state.multiplicity],
                "global_state_number": occurrence.global_state_number,
                "local_state_index": occurrence.state.local_state_index,
                "state_label": occurrence.state.state_label,
                "excitation_energy_ev": occurrence.state.energy_ev,
                "nto_filename": occurrence.nto_filename,
                "nto_path": relative_path(occurrence.nto_path, repo),
                "leading_pair_weight": leading_weight,
                "total_printed_weight": total_weight,
                "number_of_printed_pairs": len(pairs),
                "nto_character": classify_nto_character(leading_weight) if pairs else "",
                "minimum_pairs_main": main_count,
                "cumulative_weight_main": main_weight,
                "recommended_pairs_main": pair_indices(main_count),
                "minimum_pairs_si": si_count,
                "cumulative_weight_si": si_weight,
                "recommended_pairs_si": pair_indices(si_count),
                "recommendation_status": status,
                "warnings": warning_text(occurrence_warning_values(occurrence)),
                "discarded_duplicate_sources": "; ".join(occurrence.discarded_duplicate_sources),
            }
        )
    return rows


def canonical_pair_rows(repo: Path, canonical: Iterable[NTOOccurrence]) -> list[dict[str, Any]]:
    return raw_pair_rows(repo, canonical)


def cube_manifest_rows(
    repo: Path,
    canonical: Iterable[NTOOccurrence],
    cutoff: float,
) -> list[dict[str, Any]]:
    rows = []
    for occurrence in canonical:
        state = occurrence.state
        if state is None:
            continue
        pair_count, _, _ = minimum_pairs_for_cutoff(occurrence.pairs, cutoff)
        for pair in occurrence.pairs[:pair_count]:
            for role, orbital_label in [("hole", pair.hole_orbital), ("electron", pair.electron_orbital)]:
                rows.append(
                    {
                        "system": occurrence.system,
                        "system_label": SYSTEM_LABELS[occurrence.system],
                        "multiplicity": state.multiplicity,
                        "multiplicity_label": MULTIPLICITY_LABELS[state.multiplicity],
                        "state_label": state.state_label,
                        "global_state_number": occurrence.global_state_number,
                        "excitation_energy_ev": state.energy_ev,
                        "nto_filename": occurrence.nto_filename,
                        "nto_path": relative_path(occurrence.nto_path, repo),
                        "pair_index": pair.pair_index,
                        "pair_weight": pair.weight,
                        "cumulative_weight": pair.cumulative_weight,
                        "orbital_role": role,
                        "orbital_label": orbital_label,
                        "recommended_cube_basename": (
                            f"{occurrence.system}_{state.state_label}_pair{pair.pair_index}_{role}_{orbital_label}.cube"
                        ),
                    }
                )
    return rows


def parse_output_file(repo: Path, path: Path, energy_tolerance_ev: float) -> tuple[list[NTOOccurrence], dict[str, Any]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    system = infer_system(path)
    if system is None:
        raise NTOParseError(f"Could not infer PDI system for {path}")

    states = parse_state_sections(text, system, path)
    occurrences = parse_nto_blocks(text, repo, path, system, states, energy_tolerance_ev)
    normal_termination = "****ORCA TERMINATED NORMALLY****" in text
    warnings = []
    if not normal_termination:
        warnings.append("incomplete_orca_job")
    if not occurrences:
        warnings.append("no_nto_blocks")

    summary = {
        "source_output": relative_path(path, repo),
        "inferred_system": system,
        "system_label": SYSTEM_LABELS[system],
        "dedicated_multiplicity": infer_dedicated_multiplicity(path),
        "normal_termination": normal_termination,
        "number_of_singlet_states_parsed": sum(1 for state in states.values() if state.multiplicity == "singlet"),
        "number_of_triplet_states_parsed": sum(1 for state in states.values() if state.multiplicity == "triplet"),
        "number_of_nto_blocks_parsed": len(occurrences),
        "number_of_matched_nto_blocks": sum(1 for occurrence in occurrences if occurrence.state is not None),
        "missing_nto_files": sum(1 for occurrence in occurrences if not occurrence.nto_file_exists),
        "energy_mismatches": sum(1 for occurrence in occurrences if "energy_mismatch" in occurrence.warnings),
        "duplicate_conflict_status": "",
        "warnings": warning_text([*warnings, *(warning for occurrence in occurrences for warning in occurrence.warnings)]),
    }
    return occurrences, summary


def export_outputs(
    repo: Path,
    output_dir: Path,
    occurrences: list[NTOOccurrence],
    canonical: list[NTOOccurrence],
    parse_summaries: list[dict[str, Any]],
    main_cutoff: float,
    si_cutoff: float,
) -> None:
    raw_rows = raw_pair_rows(repo, occurrences)
    canonical_rows = canonical_pair_rows(repo, canonical)
    summary_rows = state_summary_rows(repo, canonical, main_cutoff, si_cutoff)
    main_manifest_rows = cube_manifest_rows(repo, canonical, main_cutoff)
    si_manifest_rows = cube_manifest_rows(repo, canonical, si_cutoff)

    write_csv(output_dir / "nto_pairs_raw.csv", raw_rows, RAW_PAIR_FIELDS)
    write_csv(output_dir / "nto_pairs.csv", canonical_rows, RAW_PAIR_FIELDS)
    write_csv(output_dir / "nto_state_summary.csv", summary_rows, STATE_SUMMARY_FIELDS)
    write_csv(output_dir / "nto_parse_summary.csv", parse_summaries, PARSE_SUMMARY_FIELDS)
    write_csv(output_dir / "nto_cube_manifest_main.csv", main_manifest_rows, CUBE_MANIFEST_FIELDS)
    write_csv(output_dir / "nto_cube_manifest_si.csv", si_manifest_rows, CUBE_MANIFEST_FIELDS)
    write_json(
        output_dir / "nto_summary.json",
        {
            "main_cutoff": main_cutoff,
            "si_cutoff": si_cutoff,
            "number_of_raw_pairs": len(raw_rows),
            "number_of_canonical_pairs": len(canonical_rows),
            "number_of_canonical_states": len(summary_rows),
            "outputs": {
                "nto_pairs_raw": "nto_pairs_raw.csv",
                "nto_pairs": "nto_pairs.csv",
                "nto_state_summary": "nto_state_summary.csv",
                "nto_parse_summary": "nto_parse_summary.csv",
                "nto_cube_manifest_main": "nto_cube_manifest_main.csv",
                "nto_cube_manifest_si": "nto_cube_manifest_si.csv",
            },
        },
    )


def update_duplicate_statuses(parse_summaries: list[dict[str, Any]], source_status: dict[str, str]) -> None:
    for summary in parse_summaries:
        status = source_status.get(summary["source_output"], "")
        if not status and int(summary["number_of_nto_blocks_parsed"]) == 0:
            status = "no_nto_blocks"
        elif not status:
            status = "no_canonical_nto"
        summary["duplicate_conflict_status"] = status


def critical_warnings(
    repo: Path,
    occurrences: Iterable[NTOOccurrence],
    parse_summaries: Iterable[dict[str, Any]],
) -> list[str]:
    critical = {"incomplete_orca_job", "unknown_state_number", "energy_mismatch", "missing_nto_file", "duplicate_conflict"}
    warnings = set()
    for occurrence in occurrences:
        for warning in occurrence.warnings:
            if warning in critical:
                warnings.add(f"{relative_path(occurrence.source_output, repo)}: {warning}")
    for summary in parse_summaries:
        for warning in str(summary.get("warnings", "")).split("; "):
            if warning in critical:
                warnings.add(f"{summary['source_output']}: {warning}")
    return sorted(warnings)


def print_console_report(
    repo: Path,
    output_dir: Path,
    canonical: list[NTOOccurrence],
    parse_summaries: list[dict[str, Any]],
) -> None:
    if not canonical:
        print("No canonical NTO states were parsed.")
    for occurrence in canonical:
        state = occurrence.state
        if state is None:
            continue
        main_count, main_weight, main_status = minimum_pairs_for_cutoff(occurrence.pairs, 0.90)
        si_count, si_weight, si_status = minimum_pairs_for_cutoff(occurrence.pairs, 0.95)
        print()
        print(f"{SYSTEM_LABELS[occurrence.system]} {state.state_label}")
        print(f"Global state: {occurrence.global_state_number}")
        print(f"Multiplicity: {MULTIPLICITY_LABELS[state.multiplicity]}")
        print(f"Energy: {state.energy_ev:.3f} eV")
        print(f"NTO file: {occurrence.nto_filename}")
        for pair in occurrence.pairs:
            print(
                f"Pair {pair.pair_index}: {pair.hole_orbital} -> {pair.electron_orbital}, "
                f"weight = {pair.weight:.8f}, cumulative = {pair.cumulative_weight:.8f}"
            )
        print(f"Main recommendation: pairs {pair_indices(main_count)} ({main_weight:.5f}, {main_status})")
        print(f"SI recommendation: pairs {pair_indices(si_count)} ({si_weight:.5f}, {si_status})")
        if occurrence.warnings:
            print(f"Warnings: {warning_text(occurrence.warnings)}")

    summary_warnings = [summary for summary in parse_summaries if summary.get("warnings")]
    if summary_warnings:
        print()
        print("Parse warnings:")
        for summary in summary_warnings:
            print(f"- {summary['source_output']}: {summary['warnings']}")

    print()
    print(f"Wrote NTO outputs under {output_dir}")


RAW_PAIR_FIELDS = [
    "system",
    "system_label",
    "source_output",
    "job_normal_termination",
    "multiplicity",
    "multiplicity_label",
    "global_state_number",
    "local_state_index",
    "state_label",
    "state_energy_hartree",
    "state_energy_ev",
    "state_energy_cm1",
    "nto_block_energy_ev",
    "nto_filename",
    "nto_path",
    "nto_file_exists",
    "nto_print_threshold",
    "pair_index",
    "hole_orbital",
    "electron_orbital",
    "weight",
    "cumulative_weight",
    "parse_warning",
]

STATE_SUMMARY_FIELDS = [
    "system",
    "system_label",
    "multiplicity",
    "multiplicity_label",
    "global_state_number",
    "local_state_index",
    "state_label",
    "excitation_energy_ev",
    "nto_filename",
    "nto_path",
    "leading_pair_weight",
    "total_printed_weight",
    "number_of_printed_pairs",
    "nto_character",
    "minimum_pairs_main",
    "cumulative_weight_main",
    "recommended_pairs_main",
    "minimum_pairs_si",
    "cumulative_weight_si",
    "recommended_pairs_si",
    "recommendation_status",
    "warnings",
    "discarded_duplicate_sources",
]

PARSE_SUMMARY_FIELDS = [
    "source_output",
    "inferred_system",
    "system_label",
    "dedicated_multiplicity",
    "normal_termination",
    "number_of_singlet_states_parsed",
    "number_of_triplet_states_parsed",
    "number_of_nto_blocks_parsed",
    "number_of_matched_nto_blocks",
    "missing_nto_files",
    "energy_mismatches",
    "duplicate_conflict_status",
    "warnings",
]

CUBE_MANIFEST_FIELDS = [
    "system",
    "system_label",
    "multiplicity",
    "multiplicity_label",
    "state_label",
    "global_state_number",
    "excitation_energy_ev",
    "nto_filename",
    "nto_path",
    "pair_index",
    "pair_weight",
    "cumulative_weight",
    "orbital_role",
    "orbital_label",
    "recommended_cube_basename",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse ORCA natural-transition-orbital calculations.")
    parser.add_argument("--project-root", "--repo", dest="project_root", type=Path, default=None)
    parser.add_argument("--calculation-root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--main-cutoff", type=float, default=0.90)
    parser.add_argument("--si-cutoff", type=float, default=0.95)
    parser.add_argument("--energy-tolerance-ev", type=float, default=0.01)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    repo = repo_root(args.project_root)
    calculation_root = resolve_cli_path(repo, args.calculation_root, Path("calculations"))
    output_dir = resolve_cli_path(repo, args.output_dir, Path("results/excited_state/nto"))

    outputs = discover_nto_outputs(calculation_root)
    if not outputs:
        raise NTOParseError(f"No NTO-associated ORCA output files found under {calculation_root}")

    occurrences: list[NTOOccurrence] = []
    parse_summaries: list[dict[str, Any]] = []
    for output in outputs:
        parsed_occurrences, summary = parse_output_file(repo, output, args.energy_tolerance_ev)
        occurrences.extend(parsed_occurrences)
        parse_summaries.append(summary)

    canonical, source_status = deduplicate_occurrences(repo, occurrences, args.energy_tolerance_ev)
    update_duplicate_statuses(parse_summaries, source_status)
    export_outputs(repo, output_dir, occurrences, canonical, parse_summaries, args.main_cutoff, args.si_cutoff)
    print_console_report(repo, output_dir, canonical, parse_summaries)

    if args.strict:
        warnings = critical_warnings(repo, occurrences, parse_summaries)
        if warnings:
            print()
            print("Strict mode failed:")
            for warning in warnings:
                print(f"- {warning}")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
