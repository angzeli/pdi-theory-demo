# 🧪 Tutorial

## 🧱 1. Build a molecular 3D structure

**Software:** ChemDraw → Avogadro

1. Draw the molecular structure in ChemDraw.
2. Save the editable source file as `.cdxml` or `.cdx` for future modification.
3. Export the structure as an MDL Molfile (`.mol`) for transfer to Avogadro.
4. Open the `.mol` file in Avogadro.
5. Check the molecular connectivity, bond orders, formal charges, stereochemistry, and hydrogen count.
6. Generate 3D coordinates if the imported structure is still two-dimensional.
7. Pre-optimize the structure using `Extensions → Optimize Geometry`.
   - Use MMFF94 when it is available and parameterized for the molecule.
   - Use UFF as a fallback for unsupported elements or fragments.
8. Inspect the pre-optimized structure carefully for:
   - unrealistic bond lengths or angles;
   - overlapping atoms;
   - incorrect protonation or formal charges;
   - distorted aromatic rings;
   - implausible conformations of flexible substituents.
9. Save an editable Avogadro copy as `.cml`.
10. Export the final pre-optimized structure as an XYZ file, for example:

```text
pdi_structure_initial.xyz
```

An XYZ file should contain:

- the total number of atoms on the first line;
- a comment on the second line;
- one element symbol and three Cartesian coordinates per atom on the remaining lines.

Example:

```text
3
Water molecule
O    0.000000    0.000000    0.000000
H    0.758602    0.000000    0.504284
H   -0.758602    0.000000    0.504284
```

A clean structure folder may look like:

```text
structures/
├── chemdraw/
│   └── pdi_structure.cdxml
├── mol/
│   └── pdi_structure.mol
├── avogadro/
│   └── pdi_structure.cml
└── xyz/
    └── pdi_structure_initial.xyz
```

Keeping each intermediate file makes the structure-generation workflow reproducible.

## ⚙️ 2. ORCA geometry optimization and frequency calculation

**Software:** ORCA

This section covers only geometry optimization and vibrational-frequency verification.

### 🗂️ 2.1 Prepare the calculation folder

Copy the XYZ file generated in Section 1 into a dedicated calculation folder, for example:

```text
calculations/
└── pdi/
    └── geometry_optimization/
        ├── pdi_structure_initial.xyz
        └── pdi_opt.inp
```

Keeping the input file and the exact XYZ structure used for the calculation in the same folder makes the job self-contained and easier to reproduce.

### 📝 2.2 Create the geometry-optimization input

Create `pdi_opt.inp` with:

```text
! r2SCAN-3c Opt TightSCF

%pal
  nprocs 4
end

%maxcore 5000

* xyzfile 0 1 pdi_structure_initial.xyz
```

The main settings are:

- `r2SCAN-3c`: an efficient composite DFT method suitable for molecular geometry optimization;
- `Opt`: requests geometry optimization;
- `TightSCF`: applies tighter SCF convergence criteria;
- `%pal`: sets the number of parallel MPI processes;
- `%maxcore`: sets the maximum memory in MB per process;
- `0 1`: specifies charge `0` and spin multiplicity `1`.

The charge and multiplicity must be changed for ions, radicals, triplet oxygen, or other open-shell systems.

### ▶️ 2.3 Run the geometry optimization

Use the validated ORCA launcher rather than invoking the ORCA binary directly:

```bash
$HOME/bin/orca611 pdi_soc.inp > pdi_opt.out 2> pdi_opt.err
```

To prevent macOS from sleeping during a long job:

```bash
caffeinate -i $HOME/bin/orca611 \
  pdi_soc.inp > pdi_opt.out 2> pdi_opt.err
```

Monitor the output in another Terminal window:

```bash
tail -f pdi_opt.out
```

Press `Ctrl+C` to stop `tail`; this does not stop ORCA.

To confirm that ORCA is still running:

```bash
ps -o pid,%cpu,%mem,etime,command -ax | grep '[o]rca'
```

### ✅ 2.4 Check the optimization result

Confirm normal termination:

```bash
grep "ORCA TERMINATED NORMALLY" pdi_opt.out
```

Expected result:

```text
****ORCA TERMINATED NORMALLY****
```

Confirm optimization convergence:

```bash
grep "THE OPTIMIZATION HAS CONVERGED" pdi_opt.out
```

Expected result:

```text
***        THE OPTIMIZATION HAS CONVERGED     ***
```

Inspect the error file:

```bash
cat pdi_opt.err
```

An empty error file is ideal.

The final optimized geometry is usually written as:

```text
pdi_opt.xyz
```

Open this file in Avogadro and check that:

- the bonding pattern is unchanged;
- no atoms overlap;
- no bonds have broken unexpectedly;
- aromatic and conjugated regions remain chemically sensible;
- flexible substituents have adopted plausible conformations.

The file `pdi_opt_trj.xyz` contains the complete optimization trajectory and can be opened in Avogadro to inspect how the structure evolved during optimization.

### 🎵 2.5 Create the frequency input

Create a separate input file named `pdi_freq.inp`:

```text
! r2SCAN-3c Freq TightSCF

%pal
  nprocs 4
end

%maxcore 5000

* xyzfile 0 1 pdi_opt.xyz
```

Use the optimized geometry from the previous step rather than the original Avogadro structure.

A frequency calculation evaluates the Hessian at the optimized geometry. It is used to determine whether the stationary point is a true local minimum.

### ▶️ 2.6 Run the frequency calculation

```bash
caffeinate -i $HOME/bin/orca611 \
  pdi_freq.inp > pdi_freq.out 2> pdi_freq.err
```

Monitor it with:

```bash
tail -f pdi_freq.out
```

Frequency calculations can remain silent for long periods while ORCA evaluates derivative integrals or solves response equations. A lack of new output does not necessarily mean that the job has stalled.

Check active ORCA processes with:

```bash
ps -o pid,%cpu,%mem,etime,command -ax | grep '[o]rca'
```

If an ORCA child process is using substantial CPU, the calculation is still active.

### ✅ 2.7 Check the frequency result

Confirm normal termination:

```bash
grep "ORCA TERMINATED NORMALLY" pdi_freq.out
```

Inspect the vibrational frequencies:

```bash
grep -A 140 "VIBRATIONAL FREQUENCIES" pdi_freq.out
```

For a true minimum:

- the first six modes should correspond to translation and rotation and should be close to zero;
- all genuine vibrational modes should be positive;
- one or two very small negative frequencies, typically only a few cm⁻¹, may arise from numerical noise or very soft motions and should be inspected rather than accepted automatically;
- a substantial imaginary frequency, such as several hundred cm⁻¹, indicates that the structure is not a local minimum.

If a significant imaginary frequency is present, visualize the corresponding mode using the ORCA Hessian file:

```bash
/Applications/Academic/orca_6_1_1/orca_pltvib \
  pdi_freq.hess 6
```

This creates an XYZ trajectory for the selected mode, which can be opened in Avogadro. Choose the mode number from the frequency table.

If the imaginary mode represents a genuine molecular distortion:

1. displace the structure slightly along the imaginary mode;
2. save the displaced structure as a new XYZ file;
3. re-optimize it with symmetry disabled;
4. repeat the frequency calculation.

A suitable re-optimization input is:

```text
! r2SCAN-3c Opt TightOpt TightSCF NoSym

%pal
  nprocs 4
end

%maxcore 5000

%geom
  Calc_Hess true
  Recalc_Hess 5
end

* xyzfile 0 1 pdi_mode_displaced.xyz
```

### 📌 2.8 Completion criteria

The Opt+Freq stage is complete only when:

1. the optimization terminates normally;
2. ORCA reports that the geometry has converged;
3. the optimized structure is chemically sensible on visual inspection;
4. the frequency calculation terminates normally;
5. no significant imaginary frequencies remain.

## 🔬 3. Single-point electronic structure calculation

**Software:** ORCA

After the optimized geometry has been verified as a true minimum by a frequency calculation, perform a higher-level single-point calculation on the validated structure. This step does not change the molecular geometry; it recalculates the electronic structure using a more suitable functional and a larger basis set.

### 🗂️ 3.1 Prepare the calculation folder

Create a dedicated folder for the electronic-structure calculation, for example:

```text
calculations/
└── pdi/
    └── electronic_structure_calculation/
        ├── pdi_opt.xyz
        └── pdi_sp.inp
```

Copy the frequency-validated optimized geometry into this folder:

```bash
cp ../geometry_optimization/pdi_opt.xyz \
  ./pdi_opt.xyz
```

Keeping a copy of the exact geometry used for the single-point calculation makes the folder self-contained and preserves the calculation provenance.

### 📝 3.2 Create the single-point input

Create `pdi_sp.inp` with:

```text
! wB97X-D4 def2-TZVP TightSCF RIJCOSX NoSym

%pal
  nprocs 8
end

%maxcore 3000

* xyzfile 0 1 pdi_opt.xyz
```

The main settings are:

- `wB97X-D4`: a range-separated hybrid density functional with D4 dispersion correction;
- `def2-TZVP`: a triple-ζ valence basis set with polarization functions;
- `TightSCF`: applies tighter SCF convergence criteria;
- `RIJCOSX`: accelerates Coulomb and exact-exchange calculations using auxiliary-basis and chain-of-spheres approximations;
- `NoSym`: disables automatic symmetry handling;
- `%pal`: sets the number of parallel MPI processes;
- `%maxcore`: sets the maximum memory in MB per process;
- `0 1`: specifies charge `0` and spin multiplicity `1`.

For a 36 GB Mac running with eight MPI processes, `%maxcore 3000` provides a practical balance between computational performance and system-memory headroom.

The charge and multiplicity must be changed for ions, radicals, triplet oxygen, or other open-shell systems.

### ▶️ 3.3 Run the calculation

Use the validated ORCA launcher:

```bash
caffeinate -i $HOME/bin/orca611 \
  pdi_sp.inp > pdi_sp.out 2> pdi_sp.err
```

Monitor the output in another Terminal window:

```bash
tail -f pdi_sp.out
```

Press `Ctrl+C` to stop `tail`; this does not stop ORCA.

To confirm that the parallel calculation is active:

```bash
ps -o pid,%cpu,%mem,etime,command -ax | grep '[o]rca'
```

Near the beginning of the output, ORCA should report the requested number of MPI processes, for example:

```text
Program running with 8 parallel MPI-processes
```

### ✅ 3.4 Check the result

Confirm normal termination:

```bash
grep "ORCA TERMINATED NORMALLY" pdi_sp.out
```

Confirm that a final electronic energy was produced:

```bash
grep "FINAL SINGLE POINT ENERGY" pdi_sp.out
```

Inspect the error file:

```bash
cat pdi_sp.err
```

An empty error file is ideal. Nonfatal MPI cleanup warnings may occasionally appear after a normally terminated calculation, but they should be recorded and distinguished from genuine calculation failures.

Search explicitly for common failure messages:

```bash
grep -Ei "fatal error|error termination|segmentation fault|MPI_ERR|aborting|library not loaded" \
  pdi_sp.out pdi_sp.err
```

The calculation should only be treated as successful if:

1. ORCA reports the expected number of MPI processes;
2. the SCF converges;
3. a final single-point energy is printed;
4. ORCA terminates normally;
5. no fatal error appears in either output stream.

### 📦 3.5 Important output files

A successful calculation typically generates:

- `pdi_sp.out`: the main human-readable ORCA output;
- `pdi_sp.err`: redirected standard error;
- `pdi_sp.gbw`: the binary wavefunction and molecular-orbital file;
- `pdi_sp.property.txt`: calculated properties, when generated;
- temporary and auxiliary files used internally by ORCA.

The `.gbw` file is the most important starting point for later orbital, electrostatic-potential, charge, and wavefunction analyses.

### 🔄 3.6 Prepare a Multiwfn-readable wavefunction file

Convert the ORCA `.gbw` file into Molden format:

```bash
/Applications/Academic/orca_6_1_1/orca_2mkl \
  pdi_sp -molden
```

This should generate:

```text
pdi_sp.molden.input
```

### 📌 3.7 Completion criteria

The single-point electronic-structure stage is complete when:

1. the calculation uses the intended frequency-validated geometry;
2. the requested parallel configuration launches correctly;
3. the SCF converges;
4. ORCA prints the final single-point energy;
5. ORCA terminates normally;
6. the `.gbw` file is successfully converted into a Multiwfn-readable format.

## 🧭 4. Multiwfn wavefunction analysis

**Software:** Multiwfn  
**Optional visualization software:** VESTA, VMD, ChimeraX, Avogadro

This section compares the ground-state electronic structures of:

- parent PDI;
- terminal-functionalized PDI.

The following analyses will be performed:

| Stage | Analysis | Available from the current wavefunctions? |
|---|---|:---:|
| 4.1 | Wavefunction validation | Yes |
| 4.2 | HOMO/LUMO analysis | Yes |
| 4.3 | Hirshfeld charges | Yes |
| 4.4 | Electrostatic potential | Yes |
| 4.5 | Mayer bond orders | Yes |
| 4.6 | QTAIM topology | Yes |
| 4.7 | ELF and LOL | Yes |
| 4.8 | Molecular DOS and fragment PDOS | Yes |

The input wavefunctions are:

```text
calculations/pdi/multiwfn_analysis/pdi_sp.molden.input
```

and:

```text
calculations/pdi_terminal_functionalized/multiwfn_analysis/pdi_terminal_functionalized_sp.molden.input
```

All analyses should use identical settings for both molecules wherever a direct comparison is intended.

### ⚠️ Important Multiwfn menu convention

At the main menu, `q` exits Multiwfn.

Many submenus expect an integer and do not accept `q`. Entering `q` in such a submenu may produce:

```text
Fortran runtime error: Bad integer for item 1 in list input
```

Use the displayed return option instead, normally:

```text
0
```

or:

```text
-10
```

Return to the main menu before entering:

```text
q
```

### 🗂️ 4.0 One-time Terminal setup

Move to the repository root:

```bash
cd "/Users/liangze/Desktop/Tsinghua 2026 Summer/pdi_h2o2_production/pdi-theory-demo"
```

Define reusable paths:

```bash
export REPO="$PWD"

export PDI="$REPO/calculations/pdi/multiwfn_analysis"
export TFPDI="$REPO/calculations/pdi_terminal_functionalized/multiwfn_analysis"

export PDI_WFN="$PDI/pdi_sp.molden.input"
export TFPDI_WFN="$TFPDI/pdi_terminal_functionalized_sp.molden.input"
```

Locate the Multiwfn executable:

```bash
export MWFN="$(command -v Multiwfn 2>/dev/null || command -v multiwfn 2>/dev/null)"
```

Check the result:

```bash
echo "$MWFN"
```

Confirm that both wavefunction files exist:

```bash
ls -lh "$PDI_WFN" "$TFPDI_WFN"
```

Create the analysis folders if necessary:

```bash
mkdir -p \
  "$PDI"/{wavefunction_validation,orbitals,charges,esp,bond_orders,qtaim,elf_lol,dos} \
  "$TFPDI"/{wavefunction_validation,orbitals,charges,esp,bond_orders,qtaim,elf_lol,dos}
```

These environment variables remain active only in the current Terminal session.

---

## ✅ 4.1 Wavefunction validation

Before generating any numerical or graphical results, confirm that Multiwfn has imported each wavefunction correctly.

### 4.1.1 Parent PDI

Enter:

```bash
cd "$PDI/wavefunction_validation"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_validation.session.log
```

After the wavefunction has loaded, inspect the startup summary.

The expected values are:

```text
Formula: H10 C24 N2 O4
Total atoms: 40
Total electrons: 200
Basis functions: 990
Occupied orbitals: 100
HOMO index: 100
LUMO index: 101
Wavefunction type: restricted closed-shell
```

At the main menu, enter:

```text
q
```

### 4.1.2 Terminal-functionalized PDI

Enter:

```bash
cd "$TFPDI/wavefunction_validation"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_validation.session.log
```

The expected values are:

```text
Total atoms: 70
Total electrons: 280
Basis functions: 1420
Occupied orbitals: 140
HOMO index: 140
LUMO index: 141
Wavefunction type: restricted closed-shell
```

At the main menu, enter:

```text
q
```

The validation stage is complete when the atom count, electron count, basis size, orbital occupation and wavefunction type are chemically consistent.

---

## 🧬 4.2 HOMO and LUMO analysis

Generate separate cube files for the HOMO and LUMO of each molecule.

Use the same:

- grid quality;
- orbital isovalue;
- viewing orientation;
- camera position;
- positive and negative orbital-lobe colours.

A visualization isovalue of approximately:

```text
±0.03 a.u.
```

is a reasonable starting point.

### 4.2.1 Parent PDI HOMO and LUMO

Enter:

```bash
cd "$PDI/orbitals"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_orbitals.session.log
```

At the Multiwfn main menu, enter:

```text
200
3
h
2
3
l
2
0
q
```

This selects `Other functions (Part 2)` → `Generate cube file for multiple orbital wavefunctions`, exports the HOMO as `h.cub`, exports the LUMO as `l.cub`, returns to the main menu, and exits.

Inspect the generated files:

```bash
ls -lt *.cub
```

Rename the generated cubes:

```bash
mv h.cub pdi_homo.cub
mv l.cub pdi_lumo.cub
```

### 4.2.2 Terminal-functionalized PDI HOMO and LUMO

Enter:

```bash
cd "$TFPDI/orbitals"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_orbitals.session.log
```

Use:

```text
200
3
h
2
3
l
2
0
q
```

Rename the generated cubes:

```bash
mv h.cub pdi_terminal_functionalized_homo.cub
mv l.cub pdi_terminal_functionalized_lumo.cub
```

### 4.2.3 Record common visualization settings

Create a record of the common visualization settings:

```bash
cat > orbital_settings.txt <<'EOF'
Orbitals: HOMO and LUMO
Grid quality: Multiwfn medium-quality grid
Recommended isovalue: +/-0.03 a.u.
Identical orientation, isovalue and rendering conventions must be used for both molecules.
EOF
```

If a scripted Multiwfn run reaches an integer menu after cube export, use `0` to return before entering `q`. Entering `q` directly in the integer prompt is the validated cause of the `Bad integer for item 1 in list input` failure.

The analysis should compare:

- HOMO and LUMO energies;
- HOMO–LUMO gap;
- orbital localization;
- contribution of the PDI core;
- contribution of the terminal groups;
- evidence for spatial separation of occupied and virtual frontier orbitals.

---

## ⚛️ 4.3 Hirshfeld charge analysis

Hirshfeld charges should be used as the primary ground-state population analysis.

Hirshfeld charges are not the same as Bader or QTAIM basin charges.

### 4.3.1 Parent PDI

Enter:

```bash
cd "$PDI/charges"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_hirshfeld.session.log
```

At the main menu, enter:

```text
7
1
1
0
q
```

This selects population analysis, Hirshfeld atomic charges, built-in sphericalized free-atom densities, returns from the population menu with `0`, and exits from the main menu with `q`.

Inspect the generated files:

```bash
ls -lh *.chg
```

The validated output is a `.chg` file containing element symbol, Cartesian coordinates and Hirshfeld charge. If your Multiwfn build uses a filename other than `charges.chg`, replace only the left-hand filename in the command below with the filename shown by `ls -lh *.chg`.

```bash
mv charges.chg pdi_hirshfeld.chg
```

If the `.chg` file was already written with the final name, leave it unchanged and verify it:

```bash
ls -lh pdi_hirshfeld.chg
```

### 4.3.2 Terminal-functionalized PDI

Enter:

```bash
cd "$TFPDI/charges"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_hirshfeld.session.log
```

Use the same settings:

```text
7
1
1
0
q
```

Rename it:

```bash
ls -lh *.chg
mv charges.chg pdi_terminal_functionalized_hirshfeld.chg
```

If your Multiwfn build uses a filename other than `charges.chg`, replace only the left-hand filename with the file shown by `ls -lh *.chg`.

For the functionalized system, the most useful comparison is not necessarily every individual atomic charge. Sum the charges over chemically meaningful fragments such as:

- PDI core;
- left terminal group;
- right terminal group;
- carbonyl oxygen atoms;
- imide nitrogen atoms.

---

## 🌈 4.4 Electrostatic-potential analysis

The electrostatic-potential stage contains three parts:

1. electron-density cube generation;
2. total ESP cube generation;
3. quantitative ESP analysis on an electron-density surface.

The full three-dimensional ESP grid is computationally expensive. On some macOS builds, this routine may use only approximately one CPU core even when several OpenMP threads are configured.

Avoid running other CPU-intensive applications during this calculation. Medium-quality ESP calculations on the larger terminal-functionalized molecule may take more than two hours on macOS, and the ESP evaluation can be effectively single-threaded.

### 4.4.0 Test cube export before a long calculation

Before running a medium-quality ESP calculation, test whether the installed Multiwfn build can successfully export an ESP cube using a low-quality grid.

Enter:

```bash
cd "$PDI/esp"
export OMP_NUM_THREADS=8
"$MWFN" "$PDI_WFN"
```

Use:

```text
5
12
1
2
0
q
```

This selects:

- grid-data analysis;
- total electrostatic potential;
- low-quality grid.

After the ESP calculation finishes, `2` exports the grid as a Gaussian cube file, `0` returns to the main menu, and `q` exits.

Check whether a file such as the following was produced:

```text
totesp.cub
```

Check:

```bash
ls -lh *.cub
```

If the low-quality export succeeds, delete the test file:

```bash
rm -f totesp.cub
```

and proceed to the medium-quality calculation.

If selecting export option `2` produces:

```text
Fortran runtime error: Bad integer for item 1 in list input
```

do not immediately repeat the multi-hour calculation with the same executable. Test the official Multiwfn macOS binary using the low-quality grid first. Only run the medium-quality ESP calculation after cube export has been validated end-to-end.

### 4.4.1 Parent PDI electron-density cube

Enter:

```bash
cd "$PDI/esp"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_density_cube.session.log
```

Use:

```text
5
1
2
2
0
q
```

This selects:

- grid-data analysis;
- electron density;
- medium-quality grid.

The second `2` exports the Gaussian cube file from the post-processing menu, `0` returns to the main menu, and `q` exits.

Rename the output:

```bash
mv density.cub pdi_density.cub
```

Confirm:

```bash
ls -lh pdi_density.cub
```

### 4.4.2 Parent PDI ESP cube

Enter:

```bash
caffeinate -i env OMP_NUM_THREADS=8 \
  "$MWFN" "$PDI_WFN"
```

Use:

```text
5
12
2
2
0
q
```

This selects:

- grid-data analysis;
- total electrostatic potential;
- medium-quality grid.

Let the calculation finish completely. The second `2` exports the ESP grid as `totesp.cub`, `0` returns to the main menu, and `q` exits.

Rename the output:

```bash
mv totesp.cub pdi_esp.cub
```

Confirm:

```bash
ls -lh pdi_density.cub pdi_esp.cub
```

The density and ESP cubes should use the same grid quality and spatial region if they will be combined or compared point-by-point.

### 4.4.3 Parent PDI quantitative surface ESP

Run:

```bash
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_esp_surface.session.log
```

At the main menu, enter:

```text
12
0
1
2
-1
-1
q
```

Start the surface analysis and wait for Multiwfn to finish evaluating the molecular surface.

The output should report quantities such as:

- global surface ESP minimum;

- global surface ESP maximum;

- local ESP extrema;

- molecular polarity index;

- positive and negative surface areas;

- polar and nonpolar surface fractions;

- average ESP and related surface statistics.

In this sequence, `12` opens quantitative molecular-surface analysis, `0` starts the analysis with the default electron-density isovalue `0.001` a.u. and mapped ESP, `1` exports `surfanalysis.txt`, `2` exports `surfanalysis.pdb`, `-1` returns to the upper surface-analysis menu, the second `-1` returns to the main menu, and `q` exits.

Back in Terminal, rename the exported files:

```bash
mv surfanalysis.txt pdi_esp_surface_statistics.txt
mv surfanalysis.pdb pdi_esp_extrema.pdb
```

Confirm that both files exist:

```bash
ls -lh \
  pdi_esp_surface_statistics.txt \
  pdi_esp_extrema.pdb
```

The text file contains the numerical surface-ESP results, while the PDB file contains the positions of the ESP extrema for visualization in programs such as VMD, VESTA, ChimeraX, or Avogadro.

### 4.4.4 Terminal-functionalized PDI density cube

Enter:

```bash
cd "$TFPDI/esp"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_density_cube.session.log
```

Use:

```text
5
1
2
2
0
q
```

Rename:

```bash
mv density.cub pdi_terminal_functionalized_density.cub
```

### 4.4.5 Terminal-functionalized PDI ESP cube

Enter:

```bash
caffeinate -i env OMP_NUM_THREADS=8 \
  "$MWFN" "$TFPDI_WFN"
```

Use:

```text
5
12
2
2
0
q
```

Rename:

```bash
mv totesp.cub pdi_terminal_functionalized_esp.cub
```

### 4.4.6 Terminal-functionalized PDI quantitative surface ESP

Run:

```bash
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_esp_surface.session.log
```

Use the same surface settings as for parent PDI:

```text
Electron-density isovalue: 0.001 a.u.
Mapped function: Total electrostatic potential
Grid spacing: 0.25 Bohr in the validated run
```

Enter:

```text
12
0
1
2
-1
-1
q
```

Start the surface analysis and wait for Multiwfn to finish evaluating the molecular surface.

The output should report quantities such as:

- global surface ESP minimum;

- global surface ESP maximum;

- local ESP extrema;

- molecular polarity index;

- positive and negative surface areas;

- polar and nonpolar surface fractions;

- average ESP and related surface statistics.

Back in Terminal, rename the exported files:

```bash
mv surfanalysis.txt pdi_terminal_functionalized_esp_surface_statistics.txt
mv surfanalysis.pdb pdi_terminal_functionalized_esp_extrema.pdb
```

Confirm that both files exist:

```bash
ls -lh \
  pdi_terminal_functionalized_esp_surface_statistics.txt \
  pdi_terminal_functionalized_esp_extrema.pdb
```

The text file contains the numerical surface-ESP results, while the PDB file contains the positions of the ESP extrema for visualization in programs such as VMD, VESTA, ChimeraX, or Avogadro.

For all final ESP figures, use the same:

- electron-density isovalue;
- colour scale;
- minimum and maximum ESP limits;
- viewing orientation;
- molecular representation.

---

## 🔗 4.5 Mayer bond-order analysis

Mayer bond orders provide a quantitative description of bonding based on the atomic-orbital overlap and density matrices.

### 4.5.1 Parent PDI

Enter:

```bash
cd "$PDI/bond_orders"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_mayer.session.log
```

At the main menu, enter:

```text
9
1
0
q
```

This selects bond-order analysis, Mayer bond-order analysis, returns to the main menu with `0`, and exits with `q`. Multiwfn writes the complete Mayer bond-order matrix as `bndmat.txt`.

Inspect:

```bash
ls -lt
```

Rename the matrix:

```bash
mv bndmat.txt pdi_mayer_bond_orders.txt
```

### 4.5.2 Terminal-functionalized PDI

Enter:

```bash
cd "$TFPDI/bond_orders"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_mayer.session.log
```

Use:

```text
9
1
0
q
```

Rename it:

```bash
mv bndmat.txt pdi_terminal_functionalized_mayer_bond_orders.txt
```

For the final comparison, extract only chemically meaningful bonds, for example:

- carbonyl C=O bonds;
- imide C–N bonds;
- selected central PDI C–C bonds;
- bonds connecting terminal groups to the PDI scaffold.

Very small off-diagonal matrix values should not automatically be interpreted as chemically meaningful bonds.

---

## 🕸️ 4.6 QTAIM topology analysis

QTAIM topology analysis identifies:

- nuclear critical points;
- bond critical points;
- ring critical points;
- cage critical points, where applicable;
- bond paths;
- electron-density properties at critical points.

This stage performs topology analysis. Full QTAIM basin integration and basin charges are a separate, more computationally demanding workflow.

### 4.6.1 Parent PDI

Enter:

```bash
cd "$PDI/qtaim"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_qtaim.session.log
```

At the main menu, enter:

```text
2
2
3
4
8
-4
4
6
0
-5
4
6
0
7
0
-10
q
```

This sequence opens topology analysis, searches critical points from nuclear positions, atom-pair midpoints and triangle centres, generates bond paths, exports `CPs.txt` and `CPs.pdb`, exports `paths.txt` and `paths.pdb`, exports all critical-point properties to `CPprop.txt`, returns to the main menu, and exits.

Rename:

```bash
mv CPs.txt pdi_CPs.txt
mv CPs.pdb pdi_CPs.pdb
mv paths.txt pdi_paths.txt
mv paths.pdb pdi_paths.pdb
mv CPprop.txt pdi_CPprop.txt
```

The critical-point property file may contain:

- electron density, `ρ`;
- Laplacian, `∇²ρ`;
- Hessian eigenvalues;
- kinetic-energy density;
- potential-energy density;
- total energy density.

Check whether the critical-point counts satisfy the topology relationship reported by Multiwfn. Failure of the topology check may indicate that one or more critical points were not located.

### 4.6.2 Terminal-functionalized PDI

Enter:

```bash
cd "$TFPDI/qtaim"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_qtaim.session.log
```

Repeat the same critical-point searches:

```text
2
2
3
4
8
-4
4
6
0
-5
4
6
0
7
0
-10
q
```

Generate bond paths and export:

```text
CPs.txt
CPs.pdb
paths.txt
paths.pdb
CPprop.txt
```

Rename:

```bash
mv CPs.txt pdi_terminal_functionalized_CPs.txt
mv CPs.pdb pdi_terminal_functionalized_CPs.pdb
mv paths.txt pdi_terminal_functionalized_paths.txt
mv paths.pdb pdi_terminal_functionalized_paths.pdb
mv CPprop.txt pdi_terminal_functionalized_CPprop.txt
```

The triangle-centre search may require more time for the 70-atom structure.

A cage-centre search is not normally necessary unless the topology check indicates missing cage critical points or the molecule contains a genuine three-dimensional cage.

---

## 💠 4.7 ELF and LOL analysis

Generate electron-localization-function and localized-orbital-locator cube files using the same medium-quality grid for both molecules.

In Multiwfn:

- real-space function `9` is ELF;
- real-space function `10` is LOL.

Both functions are dimensionless.

### 4.7.1 Parent PDI ELF

Enter:

```bash
cd "$PDI/elf_lol"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_elf.session.log
```

Use:

```text
5
9
2
2
0
q
```

This selects grid-data analysis, ELF, the medium-quality grid, exports `ELF.cub`, returns to the main menu, and exits.

Inspect:

```bash
ls -lt *.cub
```

Rename:

```bash
mv ELF.cub pdi_elf.cub
```

### 4.7.2 Parent PDI LOL

Run:

```bash
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_lol.session.log
```

Use:

```text
5
10
2
2
0
q
```

Rename:

```bash
mv LOL.cub pdi_lol.cub
```

### 4.7.3 Terminal-functionalized PDI ELF

Enter:

```bash
cd "$TFPDI/elf_lol"
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_elf.session.log
```

Use:

```text
5
9
2
2
0
q
```

Rename:

```bash
mv ELF.cub pdi_terminal_functionalized_elf.cub
```

### 4.7.4 Terminal-functionalized PDI LOL

Run:

```bash
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_lol.session.log
```

Use:

```text
5
10
2
2
0
q
```

Rename:

```bash
mv LOL.cub pdi_terminal_functionalized_lol.cub
```

The ELF and LOL analyses can illustrate:

- π-electron delocalization;
- carbonyl lone-pair localization;
- localization around terminal heteroatoms;
- changes in conjugation after terminal functionalization.

---

## 📊 4.8 Molecular DOS and fragment PDOS

For finite molecules, this analysis produces a broadened molecular-orbital density of states rather than a periodic band-structure DOS.

The initial demonstration will use element-projected fragments:

- carbon;
- nitrogen;
- oxygen;
- hydrogen.

A later, more chemically informative analysis can divide the terminal-functionalized molecule into:

- PDI core;
- terminal functional groups.

### 4.8.1 Create an atom-index helper script

Move to the repository root:

```bash
cd "$REPO"
```

Create:

```bash
mkdir -p scripts
cat > scripts/list_molden_atoms.py <<'PY'
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
PY
```

Make it executable:

```bash
chmod +x scripts/list_molden_atoms.py
```

Generate the atom-group lists:

```bash
python3 scripts/list_molden_atoms.py "$PDI_WFN" \
  | tee "$PDI/dos/atom_groups.txt"

python3 scripts/list_molden_atoms.py "$TFPDI_WFN" \
  | tee "$TFPDI/dos/atom_groups.txt"
```

Inspect:

```bash
cat "$PDI/dos/atom_groups.txt"
cat "$TFPDI/dos/atom_groups.txt"
```

Use these generated atom lists rather than manually guessing atom indices.

### 4.8.2 Parent PDI DOS and element PDOS

Enter:

```bash
cd "$PDI/dos"
"$MWFN" "$PDI_WFN" 2>&1 | tee pdi_dos.session.log
```

At the main menu, enter:

```text
10
```

In the DOS menu, configure the unit, energy range and FWHM with:

```text
8
2
-12 4 0.05
3
0.30
```

This selects eV units, an energy range of `-12` to `4` eV, a step of `0.05` eV, Gaussian broadening, and FWHM `0.30` eV.

Open the fragment-definition menu and define element fragments:

```text
-1
1
e C
q
2
e N
q
3
e O
q
4
e H
q
e
0
```

This enters fragment definition, defines fragments 1-4 using element selectors, saves each fragment with `q` inside the fragment editor, exports `DOSfrag.txt` with `e`, and returns to the DOS menu with `0`.

The element selectors should match the atom groups in:

```text
atom_groups.txt
```

For the initial demonstration, retain the default Mulliken orbital-composition method.

Generate the TDOS and fragment PDOS curves and export the numeric data:

```text
-0
3
0
-10
q
```

This draws TDOS+PDOS, exports curve and discrete-line data to text files, returns to the DOS menu, returns to the main menu, and exits.

Multiwfn 3.8 exports:

```text
DOSfrag.txt
DOS_curve.txt
DOS_line.txt
```

It does not automatically split the output into separate `*_tdos.txt` or `*_pdos_*.txt` files. `DOS_curve.txt` contains the energy column, TDOS, OPDOS, and PDOS columns for fragments 1-10. `DOS_line.txt` contains the corresponding discrete orbital-line data.

Rename the outputs descriptively:

```bash
mv DOSfrag.txt pdi_DOSfrag.txt
mv DOS_curve.txt pdi_dos_curve.txt
mv DOS_line.txt pdi_dos_line.txt
```

If the generated filenames differ between Multiwfn versions, inspect the directory first:

```bash
ls -lh
```

Record the settings:

```bash
cat > pdi_dos_settings.txt <<'EOF'
Energy unit: eV
Energy range: -12 to 4 eV
Step: 0.05 eV
Broadening: Gaussian
FWHM: 0.30 eV
Projection fragments: C, N, O and H
Orbital-composition method: Mulliken
Multiwfn numeric exports: DOS_curve.txt, DOS_line.txt, DOSfrag.txt
EOF
```

DISLIN/PDF export through the plotting menu may fail or may produce a file with a version-dependent `dislin` prefix. Treat the text export from post-processing option `3` as the reproducible result.

### 4.8.3 Terminal-functionalized PDI DOS and element PDOS

Enter:

```bash
cd "$TFPDI/dos"
cat atom_groups.txt
```

Run:

```bash
"$MWFN" "$TFPDI_WFN" 2>&1 \
  | tee pdi_terminal_functionalized_dos.session.log
```

At the main menu, enter:

```text
10
```

Use exactly the same DOS settings:

```text
8
2
-12 4 0.05
3
0.30
-1
1
e C
q
2
e N
q
3
e O
q
4
e H
q
e
0
-0
3
0
-10
q
```

Export and rename:

```bash
mv DOSfrag.txt pdi_terminal_functionalized_DOSfrag.txt
mv DOS_curve.txt pdi_terminal_functionalized_dos_curve.txt
mv DOS_line.txt pdi_terminal_functionalized_dos_line.txt
```

Record the settings:

```bash
cat > pdi_terminal_functionalized_dos_settings.txt <<'EOF'
Energy unit: eV
Energy range: -12 to 4 eV
Step: 0.05 eV
Broadening: Gaussian
FWHM: 0.30 eV
Projection fragments: C, N, O and H
Orbital-composition method: Mulliken
Multiwfn numeric exports: DOS_curve.txt, DOS_line.txt, DOSfrag.txt
EOF
```

Because the functionalized molecule contains more atoms and molecular orbitals, its absolute TDOS intensity will naturally be larger.

For final comparison figures, either:

- normalize the TDOS curves;
- compare relative fragment contributions;
- compare orbital positions rather than absolute peak heights.

---

## 📁 4.9 Expected output structure

After completing Sections 4.1–4.8, each analysis folder should contain files resembling:

```text
multiwfn_analysis/
├── wavefunction_validation/
│   └── *_validation.session.log
├── orbitals/
│   ├── *_homo.cub
│   ├── *_lumo.cub
│   └── *.session.log
├── charges/
│   ├── *_hirshfeld.chg
│   └── *_hirshfeld.session.log
├── esp/
│   ├── *_density.cub
│   ├── *_esp.cub
│   ├── *_esp_surface_statistics.txt
│   ├── *_esp_extrema.pdb
│   └── *.session.log
├── bond_orders/
│   ├── *_mayer_bond_orders.txt
│   └── *_mayer.session.log
├── qtaim/
│   ├── *_CPs.txt
│   ├── *_CPs.pdb
│   ├── *_paths.txt
│   ├── *_paths.pdb
│   ├── *_CPprop.txt
│   └── *_qtaim.session.log
├── elf_lol/
│   ├── *_elf.cub
│   ├── *_lol.cub
│   └── *.session.log
└── dos/
    ├── atom_groups.txt
    ├── *_DOSfrag.txt
    ├── *_dos_curve.txt
    ├── *_dos_line.txt
    ├── *_dos_settings.txt
    └── *_dos.session.log
```

Large regenerable files such as `.cub`, `.gbw` and `.molden.input` files may be excluded from ordinary Git tracking or managed using Git LFS.

---

## 📌 4.10 Completion criteria

The Multiwfn stage is complete when:

1. both wavefunctions load with the expected atom, electron and orbital counts;
2. HOMO and LUMO cubes have been generated for both molecules;
3. Hirshfeld charges have been exported;
4. density and ESP analyses use identical numerical settings;
5. Mayer bond-order matrices have been exported;
6. QTAIM critical points and bond paths have been generated;
7. ELF and LOL cubes have been generated;
8. TDOS and fragment PDOS data have been exported as `*_dos_curve.txt`, `*_dos_line.txt` and `*_DOSfrag.txt`;
9. all comparison plots use consistent axes, colour scales, orientations and isovalues;
10. all menu selections, software versions and analysis settings are preserved in session logs or settings files.

---

## 🧾 5. Postprocess ORCA and Multiwfn outputs

### 5.1 Run the full parser pipeline

After all ORCA and Multiwfn files have been generated, run the parser pipeline from the repository root:

```bash
cd "/Users/liangze/Desktop/Tsinghua 2026 Summer/pdi_h2o2_production/pdi-theory-demo" && /opt/anaconda3/bin/python -m scripts.postprocess.parse_all --minimum-absolute-bond-order 0.05
```

The parser runs the manifest, validation, orbital, Hirshfeld charge, cube, ESP, Mayer bond-order, QTAIM and DOS/PDOS parsers in sequence.

The expected final line is:

```text
All parsers completed successfully.
```

### 5.2 Parsed result folders

The parser writes cleaned CSV and JSON outputs under:

```text
results/ground_state/
├── manifest/
│   └── input_manifest.csv
├── validation/
│   └── wavefunction_summary.csv
├── orbitals/
│   ├── frontier_orbital_energies.csv
│   └── orbital_energies.csv
├── charges/
│   └── hirshfeld_atomic_charges.csv
├── cubes/
│   └── cube_metadata.csv
├── esp/
│   ├── esp_surface_metrics.csv
│   ├── esp_extrema_from_text.csv
│   └── esp_extrema_from_pdb.csv
├── bond_orders/
│   └── mayer_bond_orders.csv
├── qtaim/
│   ├── critical_points.csv
│   ├── critical_point_properties.csv
│   └── bond_path_points.csv
└── dos/
    ├── dos_curves.csv
    ├── dos_lines.csv
    └── dos_metadata.json
```

### 5.3 Completion criteria

The postprocessing stage is complete when:

1. the parser command ends with `All parsers completed successfully.`;
2. `frontier_orbital_energies.csv` contains populated HOMO, LUMO and gap values for both molecules;
3. `orbital_energies.csv` contains the full molecular-orbital energy table parsed from the Molden files;
4. Hirshfeld charges, ESP extrema, Mayer bond orders, QTAIM data and DOS/PDOS curves have been written under `results/ground_state/`;
5. any expected ESP warning is understood and documented rather than treated as a parser failure.

## 6. Electronic structure analysis

The electronic structure of the parent and terminal-functionalized PDI molecules is analysed using the notebooks in:

```text
analysis/
```

These notebooks use the parsed CSV files from section 5 and generate publication-ready figures under:

```text
figures/
```

Run the notebooks in numerical order, from `00_validate_parsed_data.ipynb` through `07_dos_pdos_opdos.ipynb`. The current validated analysis set contains notebooks `00` to `07`; there is no required `08_ground_state_summary.ipynb` notebook in the present workflow.

Together, these analyses compare frontier molecular orbitals, Hirshfeld atomic charges, electrostatic potential (ESP), Mayer bond orders, QTAIM topology, ELF/LOL cube statistics and DOS/PDOS/OPDOS descriptors.

### 6.1 Validation of parsed quantum-chemical data (`00_validate_parsed_data.ipynb`)

Before performing any electronic structure analysis, all parsed quantum-chemical data are automatically validated to ensure consistency and completeness.

The validation workflow verifies:

- successful parsing of ORCA output files;
- consistency between parsed orbital energies, atomic properties, and bond information;
- preservation of atom ordering across all generated datasets;
- completeness of generated CSV files required by downstream analyses.

This validation step serves as a quality-control checkpoint before higher-level analyses are performed and helps identify parsing errors or incomplete calculations at an early stage.

The notebook reads the main parsed outputs under:

```text
results/ground_state/
```

and confirms that the expected parent and terminal-functionalized PDI records are present before the remaining notebooks are used.

---

### 6.2 Frontier molecular orbital analysis (`01_frontier_orbitals.ipynb`)

Frontier molecular orbital (FMO) analysis is performed to investigate the spatial distribution and energetic characteristics of the HOMO and LUMO.

The notebook automatically:

- extracts HOMO/LUMO energies;
- calculates the HOMO–LUMO gap;
- compares orbital energies between the parent and functionalized molecules;
- exports publication-ready orbital energy tables.

The notebook reads:

```text
results/ground_state/orbitals/frontier_orbital_energies.csv
results/ground_state/orbitals/orbital_energies.csv
```

and generates:

```text
figures/orbitals/fmo_energy_level_diagram.pdf
figures/orbitals/fmo_energy_level_diagram.png
figures/orbitals/fmo_combined_figure.tiff
```

The corresponding orbital wavefunctions are visualized using the exported orbital images. For each molecule, both the HOMO and LUMO should be rendered using identical:

- isovalue;
- viewing angle;
- molecular orientation;
- rendering parameters.

This enables direct visual comparison of orbital localization before and after terminal functionalization and facilitates interpretation of changes in electron delocalization and potential charge-transfer pathways.

---

### 6.3 Hirshfeld charge analysis (`02_hirshfeld_charges.ipynb`)

Atomic charge redistribution induced by terminal functionalization is quantified using Hirshfeld population analysis.

The workflow:

- extracts Hirshfeld charges from ORCA calculations;
- maps equivalent atoms between the parent and functionalized structures;
- computes charge differences for each mapped atom;
- summarizes charge redistribution within chemically meaningful regions.

The notebook reads:

```text
results/ground_state/charges/hirshfeld_atomic_charges.csv
config/pdi_core_atom_mapping.csv
config/functionalized_atom_regions.csv
```

and generates:

```text
figures/hirshfeld_charges/matched_atom_hirshfeld_charge_comparison.pdf
figures/hirshfeld_charges/matched_atom_hirshfeld_charge_comparison.png
figures/hirshfeld_charges/matched_atom_hirshfeld_charge_shift.pdf
figures/hirshfeld_charges/matched_atom_hirshfeld_charge_shift.png
```

The resulting charge analysis reveals how electron density is redistributed between the PDI core and terminal substituents, providing quantitative insight into inductive and resonance effects introduced by functionalization.

---

### 6.4 Electrostatic potential analysis (`03_esp_analysis.ipynb`)

Electrostatic potential (ESP) analysis is performed to visualize changes in surface electrostatic characteristics resulting from terminal functionalization.

The notebook automatically:

- parses ESP cube files;
- determines consistent global colour scales;
- prepares visualization inputs;
- exports quantitative ESP statistics.

The notebook reads:

```text
results/ground_state/esp/esp_surface_metrics.csv
results/ground_state/esp/esp_extrema_from_text.csv
results/ground_state/esp/esp_extrema_from_pdb.csv
```

and writes additional summary tables:

```text
results/ground_state/esp/esp_summary.csv
results/ground_state/esp/top_10_esp_minima.csv
results/ground_state/esp/top_10_esp_maxima.csv
```

The main exported figures are:

```text
figures/esp/global_surface_esp_extrema.pdf
figures/esp/global_surface_esp_extrema.png
figures/esp/ranked_surface_esp_minima.pdf
figures/esp/ranked_surface_esp_minima.png
figures/esp/ranked_surface_esp_maxima.pdf
figures/esp/ranked_surface_esp_maxima.png
figures/esp/esp_combined_figure.tiff
```

The molecular electrostatic potential is mapped onto the electron-density isosurface using **VMD**, employing identical visualization settings for both systems, including:

- identical density isovalue;
- identical ESP colour scale;
- identical molecular orientation;
- identical rendering parameters.

The resulting ESP maps allow direct comparison of electron-rich and electron-deficient regions, enabling qualitative assessment of how terminal functionalization modifies intermolecular interaction sites and surface polarity while preserving consistency between visualizations.

---

### 6.5 Mayer bond-order analysis (`04_mayer_bond_orders.ipynb`)

Chemical bonding changes induced by terminal functionalization are investigated using Mayer bond-order analysis.

The workflow automatically:

- parses Mayer bond orders from ORCA calculations;
- maps chemically equivalent bonds between the parent and functionalized molecules;
- calculates bond-order differences;
- generates publication-ready comparison figures.

The notebook reads:

```text
results/ground_state/bond_orders/mayer_bond_orders.csv
config/selected_bonds_candidates.csv
```

and writes:

```text
results/ground_state/bond_orders/selected_mayer_bond_comparison.csv
```

The main exported figures are:

```text
figures/bond_orders/largest_selected_mayer_changes.pdf
figures/bond_orders/largest_selected_mayer_changes.png
```

Representative chemically important bonds, including carbonyl C=O bonds, imide C–N bonds, and selected aromatic C–C bonds, are compared directly between the two systems.

Both individual bond-order changes and bond-type summaries are generated, allowing identification of the bonding motifs most strongly affected by terminal functionalization. Ranking bonds by the magnitude of the Mayer bond-order change provides a concise visualization of the primary electronic perturbations introduced by functionalization while demonstrating that the conjugated PDI framework remains largely preserved.

---

### 6.6 QTAIM topology and bond-critical-point analysis (`05_qtaim_analysis.ipynb`)

QTAIM analysis is used to validate the topology of the electron density and compare bond-critical-point (BCP) descriptors between matched bonds in the parent and terminal-functionalized PDI frameworks.

The notebook reads:

```text
results/ground_state/qtaim/critical_points.csv
results/ground_state/qtaim/critical_point_properties.csv
results/ground_state/qtaim/bond_path_points.csv
config/pdi_core_atom_mapping.csv
config/selected_bonds_candidates.csv
```

The workflow:

- counts nuclear, bond and ring critical points;
- checks the Poincare-Hopf topology relation;
- summarizes BCP electron density and Laplacian distributions;
- matches BCPs through connected atom pairs rather than only by critical-point index;
- compares matched BCP descriptors for chemically selected bonds.

The notebook writes:

```text
results/ground_state/qtaim/topology_check.csv
results/ground_state/qtaim/qtaim_bcp_summary.csv
results/ground_state/qtaim/matched_bcp_comparison.csv
results/ground_state/qtaim/matched_bcp_class_summary.csv
```

and generates:

```text
figures/qtaim/critical_point_counts.pdf
figures/qtaim/critical_point_counts.png
figures/qtaim/bcp_electron_density_distribution.pdf
figures/qtaim/bcp_electron_density_distribution.png
figures/qtaim/bcp_laplacian_distribution.pdf
figures/qtaim/bcp_laplacian_distribution.png
figures/qtaim/matched_bcp_electron_density_scatter.pdf
figures/qtaim/matched_bcp_electron_density_scatter.png
figures/qtaim/matched_bcp_electron_density_change.pdf
figures/qtaim/matched_bcp_electron_density_change.png
figures/qtaim/matched_bcp_laplacian_change.pdf
figures/qtaim/matched_bcp_laplacian_change.png
```

Use the matched BCP figures to interpret local bonding changes only for bonds that are explicitly paired by the notebook. Unmatched terminal substituent bonds should be treated as functionalization-specific features rather than one-to-one perturbations of parent PDI bonds.

---

### 6.7 ELF and LOL cube analysis (`06_elf_lol_analysis.ipynb`)

ELF and LOL cube analysis compares electron-localization descriptors from the Multiwfn-generated cube files.

The notebook reads:

```text
results/ground_state/cubes/cube_metadata.csv
```

and loads the ELF/LOL cube files from the ground-state analysis folders, including:

```text
calculations/pdi/multiwfn_analysis/ground_state/elf_lol/
calculations/pdi_terminal_functionalized/multiwfn_analysis/ground_state/elf_lol/
```

The workflow:

- loads ELF and LOL volumetric grids;
- computes grid-level summary statistics;
- plots normalized voxel histograms;
- plots cumulative distribution functions;
- exports percentile tables;
- generates central z-slice visualizations for qualitative comparison.

The notebook writes:

```text
results/ground_state/cubes/elf_lol_summary.csv
results/ground_state/cubes/elf_lol_percentiles.csv
```

and generates:

```text
figures/elf_lol/pdi_elf_central_slice.pdf
figures/elf_lol/pdi_elf_central_slice.png
figures/elf_lol/pdi_lol_central_slice.pdf
figures/elf_lol/pdi_lol_central_slice.png
figures/elf_lol/pdi_terminal_functionalized_elf_central_slice.pdf
figures/elf_lol/pdi_terminal_functionalized_elf_central_slice.png
figures/elf_lol/pdi_terminal_functionalized_lol_central_slice.pdf
figures/elf_lol/pdi_terminal_functionalized_lol_central_slice.png
```

The histogram and CDF plots are quantitative comparisons of all grid voxels. The central-slice figures are qualitative two-dimensional views and should be interpreted together with the percentile summaries.

---

### 6.8 DOS, PDOS and OPDOS analysis (`07_dos_pdos_opdos.ipynb`)

DOS, PDOS and OPDOS analysis compares the electronic-state distributions of the parent and terminal-functionalized molecules.

The notebook reads:

```text
results/ground_state/dos/dos_curves.csv
results/ground_state/dos/dos_lines.csv
results/ground_state/orbitals/frontier_orbital_energies.csv
```

The workflow:

- verifies that the DOS curves and line spectra were parsed correctly;
- plots total DOS and projected DOS for each molecule;
- aligns total DOS curves relative to the HOMO energy;
- normalizes selected comparisons by atom count where appropriate;
- prepares an analysis-ready DOS table for downstream reuse.

The notebook writes:

```text
results/ground_state/dos/dos_analysis_ready.csv
```

and generates:

```text
figures/dos/pdi_tdos_pdos.pdf
figures/dos/pdi_tdos_pdos.png
figures/dos/pdi_terminal_functionalized_tdos_pdos.pdf
figures/dos/pdi_terminal_functionalized_tdos_pdos.png
figures/dos/tdos_comparison_homo_aligned_per_atom.pdf
figures/dos/tdos_comparison_homo_aligned_per_atom.png
```

The HOMO-aligned comparison is the most useful plot for comparing changes in electronic-state distribution because it removes the absolute orbital-energy offset and emphasizes changes relative to the frontier level.

---

### 6.9 Completion criteria

The electronic structure analysis stage is complete when:

1. notebooks `00_validate_parsed_data.ipynb` through `07_dos_pdos_opdos.ipynb` run without errors;
2. every expected figure folder is populated under `figures/`;
3. the additional analysis tables are present under `results/ground_state/`;
4. QTAIM matched BCP comparisons use connected-atom matching rather than critical-point order alone;
5. ELF/LOL statistics are reported from the actual cube grids rather than only from rendered images;
6. DOS comparisons include both system-specific TDOS/PDOS figures and the HOMO-aligned per-atom TDOS comparison.

---

## 7. Excited-state analysis

**Software:** ORCA, Multiwfn, Python/Jupyter

This section documents the validated excited-state workflow for parent PDI and terminal-functionalized PDI. It covers TDA singlet/triplet calculations, UV-Vis parsing, S1/T1/T2 identification, Multiwfn natural-transition-orbital analysis, Multiwfn hole-electron analysis, and ORCA spin-orbit-coupling analysis.

The validated workflow uses the `_tprint` ORCA outputs because they print excitation coefficients down to `TPrint 1e-8`, which is needed for robust Multiwfn excited-state analysis.

### 7.0 One-time Terminal setup

Start every excited-state session from the repository root and define the program paths once:

```bash
cd "/Users/liangze/Desktop/Tsinghua 2026 Summer/pdi_h2o2_production/pdi-theory-demo"

export REPO="$PWD"
export ORCA="$HOME/bin/orca611"
export ORCA_2MKL="/Applications/Academic/orca_6_1_1/orca_2mkl"
export MWFN="$(command -v Multiwfn 2>/dev/null || command -v multiwfn 2>/dev/null)"
```

Check that the executable paths resolve before running calculations:

```bash
echo "$ORCA"
echo "$ORCA_2MKL"
echo "$MWFN"
```

If `MWFN` is empty, add Multiwfn to `PATH` in the current Terminal session or call the full Multiwfn executable path directly.

Define the validated excited-state folders:

```bash
export PDI_TDA_SINGLET="$REPO/calculations/pdi/excited_state_calculations/tda_singlets"
export PDI_TDA_TRIPLET="$REPO/calculations/pdi/excited_state_calculations/tda_triplets"
export FPDI_TDA_SINGLET="$REPO/calculations/pdi_terminal_functionalized/excited_state_calculations/tda_singlets"
export FPDI_TDA_TRIPLET="$REPO/calculations/pdi_terminal_functionalized/excited_state_calculations/tda_triplets"

export PDI_NTO="$REPO/calculations/pdi/multiwfn_analysis/excited_state/nto"
export FPDI_NTO="$REPO/calculations/pdi_terminal_functionalized/multiwfn_analysis/excited_state/nto"

export PDI_HEA="$REPO/calculations/pdi/multiwfn_analysis/excited_state/hea"
export FPDI_HEA="$REPO/calculations/pdi_terminal_functionalized/multiwfn_analysis/excited_state/hea"
```

### 7.1 Excited-state electronic structure

The validated TDA calculations are:

```text
calculations/pdi/excited_state_calculations/tda_singlets/pdi_tda_singlets_tprint.inp
calculations/pdi/excited_state_calculations/tda_triplets/pdi_tda_triplets_tprint.inp
calculations/pdi_terminal_functionalized/excited_state_calculations/tda_singlets/pdi_terminal_functionalized_tda_singlets_tprint.inp
calculations/pdi_terminal_functionalized/excited_state_calculations/tda_triplets/pdi_terminal_functionalized_tda_triplets_tprint.inp
```

The shared validated method line is:

```text
! wB97X-D3 6-31+G(d,p) TightSCF RIJCOSX SMD(Water)
```

The singlet TDA block is:

```text
%tddft
  tda true
  nroots 30
  triplets false
  TPrint 1e-8
end
```

The triplet TDA block is:

```text
%tddft
  tda true
  nroots 15
  triplets true
  TPrint 1e-8
end
```

The validated resource settings were `%pal nprocs 2` for parent PDI and `%pal nprocs 3` for terminal-functionalized PDI, with `%maxcore 1200`. These settings are conservative enough to run several excited-state jobs in parallel on the validated Apple Silicon laptop.

Run each calculation from its own folder:

```bash
cd "$PDI_TDA_SINGLET"
"$ORCA" pdi_tda_singlets_tprint.inp > pdi_tda_singlets_tprint.out 2> pdi_tda_singlets_tprint.err
grep "ORCA TERMINATED NORMALLY" pdi_tda_singlets_tprint.out
"$ORCA_2MKL" pdi_tda_singlets_tprint -molden

cd "$PDI_TDA_TRIPLET"
"$ORCA" pdi_tda_triplets_tprint.inp > pdi_tda_triplets_tprint.out 2> pdi_tda_triplets_tprint.err
grep "ORCA TERMINATED NORMALLY" pdi_tda_triplets_tprint.out
"$ORCA_2MKL" pdi_tda_triplets_tprint -molden

cd "$FPDI_TDA_SINGLET"
"$ORCA" pdi_terminal_functionalized_tda_singlets_tprint.inp > pdi_terminal_functionalized_tda_singlets_tprint.out 2> pdi_terminal_functionalized_tda_singlets_tprint.err
grep "ORCA TERMINATED NORMALLY" pdi_terminal_functionalized_tda_singlets_tprint.out
"$ORCA_2MKL" pdi_terminal_functionalized_tda_singlets_tprint -molden

cd "$FPDI_TDA_TRIPLET"
"$ORCA" pdi_terminal_functionalized_tda_triplets_tprint.inp > pdi_terminal_functionalized_tda_triplets_tprint.out 2> pdi_terminal_functionalized_tda_triplets_tprint.err
grep "ORCA TERMINATED NORMALLY" pdi_terminal_functionalized_tda_triplets_tprint.out
"$ORCA_2MKL" pdi_terminal_functionalized_tda_triplets_tprint -molden
```

The parser reads the validated `_tprint.out` files and writes excited-state tables under `results/excited_state/tda_uv_vis/`:

```bash
cd "$REPO"

python -m scripts.postprocess.excited_states.parse_tda_calculations \
  --repo "$REPO" \
  --minimum-transition-weight 0.05
```

The key parsed output files are:

```text
results/excited_state/tda_uv_vis/excitation_states.csv
results/excited_state/tda_uv_vis/excitation_state_summary.csv
results/excited_state/tda_uv_vis/selected_excited_states.csv
results/excited_state/tda_uv_vis/tda_parse_summary.csv
results/excited_state/tda_uv_vis/uv_vis_interpretation_summary.csv
results/excited_state/tda_uv_vis/uv_vis_key_states.csv
```

The validated S1/T1/T2 assignments are:

| System | S1 | T1 | T2 |
|---|---:|---:|---:|
| Parent PDI | 2.763 eV | 1.617 eV | 2.949 eV |
| Terminal-functionalized PDI | 2.748 eV | 1.605 eV | 2.935 eV |

The main UV-Vis interpretation comes from S1 because it is the first bright singlet state in both molecules. T1 and T2 are the low-lying triplet states used for NTO, HEA and SOC comparisons.

### 7.2 Natural transition orbital analysis (Multiwfn)

NTO analysis is performed with Multiwfn under:

```text
calculations/pdi/multiwfn_analysis/excited_state/nto/
calculations/pdi_terminal_functionalized/multiwfn_analysis/excited_state/nto/
```

Use the `_tprint.out` files, not the older non-`_tprint` outputs. The lower `TPrint` threshold gives Multiwfn enough transition-coefficient information for the NTO decomposition.

For singlet S1, enter this Multiwfn sequence:

```text
18
6
/path/to/*_tda_singlets_tprint.out
1
3
output_S1_nto.mwfn
0
0
q
```

For triplet T1, enter:

```text
18
6
/path/to/*_tda_triplets_tprint.out
3
1
3
output_T1_nto.mwfn
0
0
q
```

For triplet T2, enter:

```text
18
6
/path/to/*_tda_triplets_tprint.out
3
2
3
output_T2_nto.mwfn
0
0
q
```

A complete set of commands is:

```bash
mkdir -p "$PDI_NTO/s1" "$PDI_NTO/t1" "$PDI_NTO/t2"
mkdir -p "$FPDI_NTO/s1" "$FPDI_NTO/t1" "$FPDI_NTO/t2"

cd "$PDI_NTO/s1"
"$MWFN" "$PDI_TDA_SINGLET/pdi_tda_singlets_tprint.molden.input"

cd "$PDI_NTO/t1"
"$MWFN" "$PDI_TDA_TRIPLET/pdi_tda_triplets_tprint.molden.input"

cd "$PDI_NTO/t2"
"$MWFN" "$PDI_TDA_TRIPLET/pdi_tda_triplets_tprint.molden.input"

cd "$FPDI_NTO/s1"
"$MWFN" "$FPDI_TDA_SINGLET/pdi_terminal_functionalized_tda_singlets_tprint.molden.input"

cd "$FPDI_NTO/t1"
"$MWFN" "$FPDI_TDA_TRIPLET/pdi_terminal_functionalized_tda_triplets_tprint.molden.input"

cd "$FPDI_NTO/t2"
"$MWFN" "$FPDI_TDA_TRIPLET/pdi_terminal_functionalized_tda_triplets_tprint.molden.input"
```

When Multiwfn asks for the ORCA output path, use the absolute path to the matching `_tprint.out` file. Example for parent PDI S1:

```text
/Users/liangze/Desktop/Tsinghua 2026 Summer/pdi_h2o2_production/pdi-theory-demo/calculations/pdi/excited_state_calculations/tda_singlets/pdi_tda_singlets_tprint.out
```

After each `.mwfn` file is generated, export the NTO hole and electron cubes from that `.mwfn` file. Use `100 100 100` as the grid resolution. The `150 150 150` grid generated cube files that were too large for practical VMD rendering on the validated laptop.

For each NTO orbital index, enter:

```text
5
4
ORBITAL_INDEX
4
100 100 100
2
0
q
```

Multiwfn writes the cube as:

```text
MOvalue.cub
```

Rename it immediately before exporting the next orbital:

```bash
mv MOvalue.cub pdi_T2_nto_pair2_electron.cub
```

Use this validated NTO cube index table:

| System | State | Pair | Hole index | Electron index |
|---|---|---:|---:|---:|
| Parent PDI | S1 | 1 | 100 | 101 |
| Parent PDI | T1 | 1 | 100 | 101 |
| Parent PDI | T2 | 1 | 100 | 101 |
| Parent PDI | T2 | 2 | 99 | 102 |
| Functionalized PDI | S1 | 1 | 140 | 141 |
| Functionalized PDI | T1 | 1 | 140 | 141 |
| Functionalized PDI | T2 | 1 | 140 | 141 |
| Functionalized PDI | T2 | 2 | 139 | 142 |

The validated leading NTO weights are:

| System | State | Pair 1 weight | Pair 2 weight | Main rendering |
|---|---|---:|---:|---|
| Parent PDI | S1 | 0.95548 | - | Pair 1 |
| Parent PDI | T1 | 0.90821 | - | Pair 1 |
| Parent PDI | T2 | 0.49697 | 0.43201 | Pairs 1 and 2 |
| Terminal-functionalized PDI | S1 | 0.95515 | - | Pair 1 |
| Terminal-functionalized PDI | T1 | 0.90868 | - | Pair 1 |
| Terminal-functionalized PDI | T2 | 0.49176 | 0.43827 | Pairs 1 and 2 |

The practical rendering rule used here is:

- if pair 1 weight is at least `0.90`, render pair 1 only;
- if pair 1 is between `0.75` and `0.90`, add pairs until cumulative weight is at least `0.90`;
- if a secondary pair is at least `0.10`, render it;
- if a secondary pair is between `0.05` and `0.10`, report it numerically or render it in supporting information when spatially distinctive;
- if a secondary pair is below `0.05`, usually report it numerically but do not render it.

After cube generation, parse the NTO analysis:

```bash
cd "$REPO"

python -m scripts.postprocess.excited_states.parse_nto_calculations \
  --project-root "$REPO" \
  --calculation-root calculations \
  --output-dir results/excited_state/nto_multiwfn \
  --main-cutoff 0.90 \
  --si-cutoff 0.95 \
  --strict
```

The main parsed outputs are:

```text
results/excited_state/nto_multiwfn/nto_cube_manifest_main.csv
results/excited_state/nto_multiwfn/nto_cube_manifest_si.csv
results/excited_state/nto_multiwfn/nto_pairs.csv
results/excited_state/nto_multiwfn/nto_state_summary.csv
results/excited_state/nto_multiwfn/nto_summary.json
```

### 7.3 Hole-electron analysis (HEA)

HEA is also performed in Multiwfn and is stored under:

```text
calculations/pdi/multiwfn_analysis/excited_state/hea/
calculations/pdi_terminal_functionalized/multiwfn_analysis/excited_state/hea/
```

Use excitation analysis option `1`, not NTO option `6`.

For singlet S1, enter:

```text
18
1
/path/to/*_tda_singlets_tprint.out
1
2
10
1
11
1
0
0
q
```

For triplet T1 or T2, enter:

```text
18
1
/path/to/*_tda_triplets_tprint.out
3
STATE_INDEX
1
2
10
1
11
1
0
0
q
```

Use:

```text
STATE_INDEX = 1 for T1
STATE_INDEX = 2 for T2
```

A complete set of folders and launch commands is:

```bash
mkdir -p "$PDI_HEA/s1" "$PDI_HEA/t1" "$PDI_HEA/t2"
mkdir -p "$FPDI_HEA/s1" "$FPDI_HEA/t1" "$FPDI_HEA/t2"

cd "$PDI_HEA/s1"
"$MWFN" "$PDI_TDA_SINGLET/pdi_tda_singlets_tprint.molden.input"

cd "$PDI_HEA/t1"
"$MWFN" "$PDI_TDA_TRIPLET/pdi_tda_triplets_tprint.molden.input"

cd "$PDI_HEA/t2"
"$MWFN" "$PDI_TDA_TRIPLET/pdi_tda_triplets_tprint.molden.input"

cd "$FPDI_HEA/s1"
"$MWFN" "$FPDI_TDA_SINGLET/pdi_terminal_functionalized_tda_singlets_tprint.molden.input"

cd "$FPDI_HEA/t1"
"$MWFN" "$FPDI_TDA_TRIPLET/pdi_terminal_functionalized_tda_triplets_tprint.molden.input"

cd "$FPDI_HEA/t2"
"$MWFN" "$FPDI_TDA_TRIPLET/pdi_terminal_functionalized_tda_triplets_tprint.molden.input"
```

Multiwfn writes:

```text
hole.cub
electron.cub
```

Rename these immediately after each state. Examples:

```bash
mv hole.cub pdi_T1_hole.cub
mv electron.cub pdi_T1_electron.cub
```

For terminal-functionalized PDI, use the full system prefix:

```bash
mv hole.cub pdi_terminal_functionalized_T1_hole.cub
mv electron.cub pdi_terminal_functionalized_T1_electron.cub
```

If Multiwfn reports that configurations with absolute coefficient below `0.01000` are ignored, that is the normal excitation-analysis coefficient threshold and not an error. The `_tprint` ORCA files are still required because they preserve enough small transition terms for reproducible postprocessing.

### 7.4 Spin-orbit coupling analysis

SOC calculations are run directly in ORCA from:

```text
calculations/pdi/excited_state_calculations/soc/pdi_soc.inp
calculations/pdi_terminal_functionalized/excited_state_calculations/soc/pdi_terminal_functionalized_soc.inp
```

The validated SOC input core is:

```text
%tddft
    TDA true
    NRoots 10
    DoSOC true
    TPrint 1e-8
end
```

Run and check parent PDI:

```bash
cd "$REPO/calculations/pdi/excited_state_calculations/soc"
"$ORCA" pdi_soc.inp > pdi_soc.out 2> pdi_soc.err
grep "ORCA TERMINATED NORMALLY" pdi_soc.out
```

Run and check terminal-functionalized PDI:

```bash
cd "$REPO/calculations/pdi_terminal_functionalized/excited_state_calculations/soc"
"$ORCA" pdi_terminal_functionalized_soc.inp > pdi_terminal_functionalized_soc.out 2> pdi_terminal_functionalized_soc.err
grep "ORCA TERMINATED NORMALLY" pdi_terminal_functionalized_soc.out
```

Parse the SOC outputs from the repository root. Do not leave trailing spaces after the backslashes:

```bash
cd "$REPO"

python scripts/postprocess/excited_states/parse_soc.py \
  calculations/pdi/excited_state_calculations/soc/pdi_soc.out \
  --singlet 1 \
  --csv results/excited_state/soc/pdi_soc.csv

python scripts/postprocess/excited_states/parse_soc.py \
  calculations/pdi_terminal_functionalized/excited_state_calculations/soc/pdi_terminal_functionalized_soc.out \
  --singlet 1 \
  --csv results/excited_state/soc/pdi_terminal_functionalized_soc.csv
```

The validated SOC analysis focuses on:

- SOC matrix heatmaps;
- S1-triplet SOC bar charts;
- S1-triplet SOC versus energy-gap screening;
- the screening descriptor `SOC^2 / DeltaE^2`.

The validated low-lying S1-triplet SOC values are small:

| System | S1-T1 SOC | S1-T2 SOC | Strongest S1-Tn channel |
|---|---:|---:|---|
| Parent PDI | 0.090 cm^-1 | 0.014 cm^-1 | S1-T9, 10.63 cm^-1 |
| Terminal-functionalized PDI | 0.052 cm^-1 | 0.050 cm^-1 | S1-T10, 10.58 cm^-1 |

### 7.5 Excited-state notebooks and figures

The validated excited-state notebooks are:

```text
analysis/08_excitation_state_uv_vis_analysis.ipynb
analysis/09_natural_transition_orbital_analysis.ipynb
analysis/10_spin_orbit_coupling_analysis.ipynb
```

Optional notebook execution commands:

```bash
cd "$REPO"

jupyter nbconvert --to notebook --execute \
  analysis/08_excitation_state_uv_vis_analysis.ipynb \
  --output /tmp/08_excited_state_uv_vis.executed.ipynb \
  --ExecutePreprocessor.timeout=600

jupyter nbconvert --to notebook --execute \
  analysis/09_natural_transition_orbital_analysis.ipynb \
  --output /tmp/09_nto.executed.ipynb \
  --ExecutePreprocessor.timeout=600

jupyter nbconvert --to notebook --execute \
  analysis/10_spin_orbit_coupling_analysis.ipynb \
  --output /tmp/10_soc.executed.ipynb \
  --ExecutePreprocessor.timeout=600
```

The generated figure folders are:

```text
figures/excited_state_uv_vis/
figures/nto/
figures/hea/
figures/soc/
```

The validated UV-Vis and excited-state figures include:

```text
figures/excited_state_uv_vis/simulated_uv_vis_comparison_normalized.pdf
figures/excited_state_uv_vis/pdi_singlet_triplet_energy_levels.pdf
figures/excited_state_uv_vis/pdi_terminal_functionalized_singlet_triplet_energy_levels.pdf
```

The validated NTO and HEA figures include:

```text
figures/nto/pdi/pdi_nto_combined_figure.tiff
figures/nto/pdi_terminal_functionalized/pdi_terminal_functionalized_nto_combined_figure.tiff
figures/hea/hea_combined_figure.tiff
```

The validated SOC figures include:

```text
figures/soc/pdi_soc_matrix_heatmap.pdf
figures/soc/pdi_terminal_functionalized_soc_matrix_heatmap.pdf
figures/soc/s1_triplet_soc_comparison.pdf
figures/soc/s1_triplet_soc_vs_energy_gap.pdf
figures/soc/dominant_s1_soc_channels.pdf
```

Do not describe a SOC-arrow energy-level figure as a validated output unless it is generated and added later.

### 7.6 Expected output structure

After the excited-state workflow is complete, the relevant output tree should include:

```text
calculations/
├── pdi/
│   ├── excited_state_calculations/
│   │   ├── tda_singlets/
│   │   ├── tda_triplets/
│   │   └── soc/
│   └── multiwfn_analysis/
│       └── excited_state/
│           ├── hea/
│           │   ├── s1/
│           │   ├── t1/
│           │   └── t2/
│           └── nto/
│               ├── s1/
│               ├── t1/
│               └── t2/
└── pdi_terminal_functionalized/
    ├── excited_state_calculations/
    │   ├── tda_singlets/
    │   ├── tda_triplets/
    │   └── soc/
    └── multiwfn_analysis/
        └── excited_state/
            ├── hea/
            │   ├── s1/
            │   ├── t1/
            │   └── t2/
            └── nto/
                ├── s1/
                ├── t1/
                └── t2/

results/
└── excited_state/
    ├── nto_multiwfn/
    ├── soc/
    └── tda_uv_vis/

figures/
├── excited_state_uv_vis/
├── hea/
├── nto/
└── soc/
```

The key analysis-ready CSV files are:

```text
results/excited_state/tda_uv_vis/selected_excited_states.csv
results/excited_state/tda_uv_vis/excitation_state_summary.csv
results/excited_state/nto_multiwfn/nto_state_summary.csv
results/excited_state/nto_multiwfn/nto_cube_manifest_main.csv
results/excited_state/soc/s1_triplet_soc_summary.csv
results/excited_state/soc/s1_triplet_soc_energy_gap_screening.csv
```

### 7.7 Completion criteria

The excited-state analysis stage is complete when:

1. all four `_tprint` TDA calculations terminate normally;
2. all four `_tprint` calculations have corresponding `.molden.input` files from `orca_2mkl`;
3. `parse_tda_calculations.py` completes and writes `results/excited_state/tda_uv_vis/selected_excited_states.csv`;
4. S1, T1 and T2 are identified for both parent PDI and terminal-functionalized PDI;
5. Multiwfn NTO `.mwfn` files and cube files are generated for S1, T1 and T2;
6. T2 includes the two dominant NTO pairs for both molecules;
7. HEA hole and electron cubes are generated for S1, T1 and T2;
8. SOC calculations terminate normally and `parse_soc.py` writes both system-level SOC CSV files;
9. notebooks `08_excitation_state_uv_vis_analysis.ipynb`, `09_natural_transition_orbital_analysis.ipynb` and `10_spin_orbit_coupling_analysis.ipynb` run without errors;
10. the expected excited-state figure folders are populated under `figures/`.
