from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.postprocess.excited_states import parse_nto_calculations as nto

NORMAL_TERMINATION = "****ORCA TERMINATED NORMALLY****"


def project_root() -> Path:
    root = Path(tempfile.mkdtemp())
    (root / "pyproject.toml").write_text("[project]\nname='fixture'\n", encoding="utf-8")
    return root


def pdi_output_path(root: Path, folder: str, filename: str) -> Path:
    path = root / "calculations" / "pdi" / "excited_state_calculations" / folder / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def state_sections() -> str:
    return """
TD-DFT/TDA EXCITED STATES (SINGLETS)

STATE  1:  E=   0.101556 au      2.763 eV    22289.0 cm**-1 <S**2> =   0.000000 Mult 1
STATE  2:  E=   0.141875 au      3.861 eV    31138.0 cm**-1 <S**2> =   0.000000 Mult 1

TD-DFT/TDA EXCITED STATES (TRIPLETS)

STATE 31:  E=   0.059425 au      1.617 eV    13042.2 cm**-1 <S**2> =   2.000000 Mult 3
STATE 32:  E=   0.108380 au      2.949 eV    23786.6 cm**-1 <S**2> =   2.000000 Mult 3
"""


def nto_block(state_number: int, filename: str, energy_ev: float = 1.617, pairs: str | None = None) -> str:
    if pairs is None:
        pairs = """
    99a -> 100a  : n=  0.82461000
    98a -> 101a  : n=  0.09563000
    97a -> 102a  : n=  0.05213000
"""
    return f"""
------------------------------------------
NATURAL TRANSITION ORBITALS FOR STATE {state_number:4d}
------------------------------------------

Natural Transition Orbitals were saved in {filename}
Threshold for printing occupation numbers 0.000100

 E=   0.059425 au      {energy_ev:.3f} eV    13042.2 cm**-1
{pairs}
"""


class NTOParserTests(unittest.TestCase):
    def test_singlet_state_one_maps_to_s1(self) -> None:
        states = nto.parse_state_sections(state_sections(), "pdi", Path("pdi_nto_singlet.out"))

        self.assertEqual(states[1].multiplicity, "singlet")
        self.assertEqual(states[1].state_label, "S1")

    def test_triplet_global_state_31_maps_to_t1(self) -> None:
        states = nto.parse_state_sections(state_sections(), "pdi", Path("pdi_nto_triplet.out"))

        self.assertEqual(states[31].multiplicity, "triplet")
        self.assertEqual(states[31].state_label, "T1")
        self.assertAlmostEqual(states[31].energy_ev, 1.617)

    def test_s31_nto_filename_is_not_treated_as_singlet_31(self) -> None:
        root = project_root()
        path = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        text = state_sections() + nto_block(31, "pdi_nto_triplet.s31.nto") + NORMAL_TERMINATION
        path.write_text(text, encoding="utf-8")
        (path.parent / "pdi_nto_triplet.s31.nto").write_text("fixture", encoding="utf-8")

        occurrences, _ = nto.parse_output_file(root, path, 0.01)

        self.assertEqual(occurrences[0].state.state_label, "T1")
        self.assertEqual(occurrences[0].nto_filename, "pdi_nto_triplet.s31.nto")

    def test_exact_nto_filename_extraction_and_energy_match(self) -> None:
        root = project_root()
        path = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        text = state_sections() + nto_block(31, "custom_name.s31.nto") + NORMAL_TERMINATION
        path.write_text(text, encoding="utf-8")
        (path.parent / "custom_name.s31.nto").write_text("fixture", encoding="utf-8")

        occurrences, _ = nto.parse_output_file(root, path, 0.01)

        self.assertEqual(occurrences[0].nto_filename, "custom_name.s31.nto")
        self.assertNotIn("energy_mismatch", occurrences[0].warnings)

    def test_energy_mismatch_warning(self) -> None:
        root = project_root()
        path = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        text = state_sections() + nto_block(31, "pdi_nto_triplet.s31.nto", energy_ev=1.700) + NORMAL_TERMINATION
        path.write_text(text, encoding="utf-8")
        (path.parent / "pdi_nto_triplet.s31.nto").write_text("fixture", encoding="utf-8")

        occurrences, _ = nto.parse_output_file(root, path, 0.01)

        self.assertIn("energy_mismatch", occurrences[0].warnings)

    def test_cumulative_weight_recommendations(self) -> None:
        pairs = [
            nto.NTOPair(1, "99a", "100a", 0.82461, 0.82461),
            nto.NTOPair(2, "98a", "101a", 0.09563, 0.92024),
            nto.NTOPair(3, "97a", "102a", 0.05213, 0.97237),
        ]

        self.assertEqual(nto.minimum_pairs_for_cutoff(pairs, 0.90), (2, 0.92024, "reached"))
        self.assertEqual(nto.minimum_pairs_for_cutoff(pairs, 0.95), (3, 0.97237, "reached"))

    def test_total_printed_weight_below_cutoff(self) -> None:
        pairs = [
            nto.NTOPair(1, "99a", "100a", 0.40, 0.40),
            nto.NTOPair(2, "98a", "101a", 0.20, 0.60),
        ]

        self.assertEqual(nto.minimum_pairs_for_cutoff(pairs, 0.90), (2, 0.60, "printed_weights_below_cutoff"))

    def test_missing_nto_file_warning(self) -> None:
        root = project_root()
        path = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        path.write_text(state_sections() + nto_block(31, "missing.s31.nto") + NORMAL_TERMINATION, encoding="utf-8")

        occurrences, _ = nto.parse_output_file(root, path, 0.01)

        self.assertIn("missing_nto_file", occurrences[0].warnings)

    def test_duplicate_deduplication_prefers_dedicated_output(self) -> None:
        root = project_root()
        singlet = pdi_output_path(root, "nto_singlet", "pdi_nto_singlet.out")
        triplet_named = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        pair_text = "    99a -> 100a  : n=  0.95548000\n"
        singlet.write_text(state_sections() + nto_block(1, "pdi_nto_singlet.s1.nto", 2.763, pair_text) + NORMAL_TERMINATION, encoding="utf-8")
        triplet_named.write_text(
            state_sections() + nto_block(1, "pdi_nto_triplet.s1.nto", 2.763, pair_text) + NORMAL_TERMINATION,
            encoding="utf-8",
        )
        (singlet.parent / "pdi_nto_singlet.s1.nto").write_text("fixture", encoding="utf-8")
        (triplet_named.parent / "pdi_nto_triplet.s1.nto").write_text("fixture", encoding="utf-8")

        occurrences = []
        for path in [singlet, triplet_named]:
            parsed, _ = nto.parse_output_file(root, path, 0.01)
            occurrences.extend(parsed)
        canonical, _ = nto.deduplicate_occurrences(root, occurrences, 0.01)

        self.assertEqual(len(canonical), 1)
        self.assertEqual(canonical[0].nto_filename, "pdi_nto_singlet.s1.nto")
        self.assertEqual(canonical[0].duplicate_status, "duplicate_selected")

    def test_conflicting_duplicate_detection(self) -> None:
        root = project_root()
        first = pdi_output_path(root, "nto_singlet", "pdi_nto_singlet.out")
        second = pdi_output_path(root, "nto_copy", "pdi_nto_copy.out")
        first.write_text(
            state_sections() + nto_block(1, "first.s1.nto", 2.763, "    99a -> 100a  : n=  0.95548000\n") + NORMAL_TERMINATION,
            encoding="utf-8",
        )
        second.write_text(
            state_sections() + nto_block(1, "second.s1.nto", 2.763, "    98a -> 101a  : n=  0.95548000\n") + NORMAL_TERMINATION,
            encoding="utf-8",
        )
        (first.parent / "first.s1.nto").write_text("fixture", encoding="utf-8")
        (second.parent / "second.s1.nto").write_text("fixture", encoding="utf-8")

        occurrences = []
        for path in [first, second]:
            parsed, _ = nto.parse_output_file(root, path, 0.01)
            occurrences.extend(parsed)
        canonical, _ = nto.deduplicate_occurrences(root, occurrences, 0.01)

        self.assertIn("duplicate_conflict", canonical[0].warnings)

    def test_incomplete_orca_job_handling(self) -> None:
        root = project_root()
        path = pdi_output_path(root, "nto_triplet", "pdi_nto_triplet.out")
        path.write_text(state_sections() + nto_block(31, "pdi_nto_triplet.s31.nto"), encoding="utf-8")
        (path.parent / "pdi_nto_triplet.s31.nto").write_text("fixture", encoding="utf-8")

        occurrences, summary = nto.parse_output_file(root, path, 0.01)

        self.assertFalse(summary["normal_termination"])
        self.assertIn("incomplete_orca_job", occurrences[0].warnings)


if __name__ == "__main__":
    unittest.main()
