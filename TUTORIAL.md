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

    `pdi_structure_initial.xyz`

An XYZ file should contain:

- the total number of atoms on the first line;
- a comment on the second line;
- one element symbol and three Cartesian coordinates per atom on the remaining lines.

Example:

    3
    Water molecule
    O    0.000000    0.000000    0.000000
    H    0.758602    0.000000    0.504284
    H   -0.758602    0.000000    0.504284

A clean structure folder may look like:

    structures/
    ├── chemdraw/
    │   └── pdi_structure.cdxml
    ├── mol/
    │   └── pdi_structure.mol
    ├── avogadro/
    │   └── pdi_structure.cml
    └── xyz/
        └── pdi_structure_initial.xyz

Keeping each intermediate file makes the structure-generation workflow reproducible.

## ⚙️ 2. ORCA geometry optimization and frequency calculation

**Software:** ORCA

This section covers only geometry optimization and vibrational-frequency verification.

### 🗂️ 2.1 Prepare the calculation folder

Copy the XYZ file generated in Section 1 into a dedicated calculation folder, for example:

    calculations/
    └── pdi/
        └── geometry_optimization/
            ├── pdi_structure_initial.xyz
            └── pdi_opt.inp

Keeping the input file and the exact XYZ structure used for the calculation in the same folder makes the job self-contained and easier to reproduce.

### 📝 2.2 Create the geometry-optimization input

Create `pdi_opt.inp` with:

    ! r2SCAN-3c Opt TightSCF

    %pal
      nprocs 4
    end

    %maxcore 5000

    * xyzfile 0 1 pdi_structure_initial.xyz

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

    $HOME/bin/orca611 pdi_opt.inp > pdi_opt.out 2> pdi_opt.err

To prevent macOS from sleeping during a long job:

    caffeinate -i $HOME/bin/orca611 \
      pdi_opt.inp > pdi_opt.out 2> pdi_opt.err

Monitor the output in another Terminal window:

    tail -f pdi_opt.out

Press `Ctrl+C` to stop `tail`; this does not stop ORCA.

To confirm that ORCA is still running:

    ps -o pid,%cpu,%mem,etime,command -ax | grep '[o]rca'

### ✅ 2.4 Check the optimization result

Confirm normal termination:

    grep "ORCA TERMINATED NORMALLY" pdi_opt.out

Expected result:

    ****ORCA TERMINATED NORMALLY****

Confirm optimization convergence:

    grep "THE OPTIMIZATION HAS CONVERGED" pdi_opt.out

Expected result:

    ***        THE OPTIMIZATION HAS CONVERGED     ***

Inspect the error file:

    cat pdi_opt.err

An empty error file is ideal.

The final optimized geometry is usually written as:

    pdi_opt.xyz

Open this file in Avogadro and check that:

- the bonding pattern is unchanged;
- no atoms overlap;
- no bonds have broken unexpectedly;
- aromatic and conjugated regions remain chemically sensible;
- flexible substituents have adopted plausible conformations.

The file `pdi_opt_trj.xyz` contains the complete optimization trajectory and can be opened in Avogadro to inspect how the structure evolved during optimization.

### 🎵 2.5 Create the frequency input

Create a separate input file named `pdi_freq.inp`:

    ! r2SCAN-3c Freq TightSCF

    %pal
      nprocs 4
    end

    %maxcore 5000

    * xyzfile 0 1 pdi_opt.xyz

Use the optimized geometry from the previous step rather than the original Avogadro structure.

A frequency calculation evaluates the Hessian at the optimized geometry. It is used to determine whether the stationary point is a true local minimum.

### ▶️ 2.6 Run the frequency calculation

    caffeinate -i $HOME/bin/orca611 \
      pdi_freq.inp > pdi_freq.out 2> pdi_freq.err

Monitor it with:

    tail -f pdi_freq.out

Frequency calculations can remain silent for long periods while ORCA evaluates derivative integrals or solves response equations. A lack of new output does not necessarily mean that the job has stalled.

Check active ORCA processes with:

    ps -o pid,%cpu,%mem,etime,command -ax | grep '[o]rca'

If an ORCA child process is using substantial CPU, the calculation is still active.

### ✅ 2.7 Check the frequency result

Confirm normal termination:

    grep "ORCA TERMINATED NORMALLY" pdi_freq.out

Inspect the vibrational frequencies:

    grep -A 140 "VIBRATIONAL FREQUENCIES" pdi_freq.out

For a true minimum:

- the first six modes should correspond to translation and rotation and should be close to zero;
- all genuine vibrational modes should be positive;
- one or two very small negative frequencies, typically only a few cm⁻¹, may arise from numerical noise or very soft motions and should be inspected rather than accepted automatically;
- a substantial imaginary frequency, such as several hundred cm⁻¹, indicates that the structure is not a local minimum.

If a significant imaginary frequency is present, visualize the corresponding mode using the ORCA Hessian file:

    /Applications/Academic/orca_6_1_1/orca_pltvib \
      pdi_freq.hess 6

This creates an XYZ trajectory for the selected mode, which can be opened in Avogadro. Choose the mode number from the frequency table.

If the imaginary mode represents a genuine molecular distortion:

1. displace the structure slightly along the imaginary mode;
2. save the displaced structure as a new XYZ file;
3. re-optimize it with symmetry disabled;
4. repeat the frequency calculation.

A suitable re-optimization input is:

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

```text
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

```bash
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

    calculations/pdi/multiwfn_analysis/pdi_sp.molden.input

and:

    calculations/pdi_terminal_functionalized/multiwfn_analysis/
    pdi_terminal_functionalized_sp.molden.input

All analyses should use identical settings for both molecules wherever a direct comparison is intended.

### ⚠️ Important Multiwfn menu convention

At the main menu, `q` exits Multiwfn.

Many submenus expect an integer and do not accept `q`. Entering `q` in such a submenu may produce:

    Fortran runtime error: Bad integer for item 1 in list input

Use the displayed return option instead, normally:

    0

or:

    -10

Return to the main menu before entering:

    q

### 🗂️ 4.0 One-time Terminal setup

Move to the repository root:

    cd "/Users/liangze/Desktop/Tsinghua 2026 Summer/pdi_h2o2_production/pdi-theory-demo"

Define reusable paths:

    export REPO="$PWD"

    export PDI="$REPO/calculations/pdi/multiwfn_analysis"
    export TFPDI="$REPO/calculations/pdi_terminal_functionalized/multiwfn_analysis"

    export PDI_WFN="$PDI/pdi_sp.molden.input"
    export TFPDI_WFN="$TFPDI/pdi_terminal_functionalized_sp.molden.input"

Locate the Multiwfn executable:

    export MWFN="$(command -v Multiwfn 2>/dev/null || command -v multiwfn 2>/dev/null)"

Check the result:

    echo "$MWFN"

Confirm that both wavefunction files exist:

    ls -lh "$PDI_WFN" "$TFPDI_WFN"

Create the analysis folders if necessary:

    mkdir -p \
      "$PDI"/{wavefunction_validation,orbitals,charges,esp,bond_orders,qtaim,elf_lol,dos,cubes} \
      "$TFPDI"/{wavefunction_validation,orbitals,charges,esp,bond_orders,qtaim,elf_lol,dos,cubes}

These environment variables remain active only in the current Terminal session.

---

## ✅ 4.1 Wavefunction validation

Before generating any numerical or graphical results, confirm that Multiwfn has imported each wavefunction correctly.

### 4.1.1 Parent PDI

Enter:

    cd "$PDI/wavefunction_validation"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_validation.session.log

After the wavefunction has loaded, inspect the startup summary.

The expected values are:

    Formula: H10 C24 N2 O4
    Total atoms: 40
    Total electrons: 200
    Basis functions: 990
    Occupied orbitals: 100
    HOMO index: 100
    LUMO index: 101
    Wavefunction type: restricted closed-shell

At the main menu, enter:

    q

### 4.1.2 Terminal-functionalized PDI

Enter:

    cd "$TFPDI/wavefunction_validation"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_validation.session.log

The expected values are:

    Total atoms: 70
    Total electrons: 280
    Basis functions: 1420
    Occupied orbitals: 140
    HOMO index: 140
    LUMO index: 141
    Wavefunction type: restricted closed-shell

At the main menu, enter:

    q

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

    ±0.03 a.u.

is a reasonable starting point.

### 4.2.1 Parent PDI HOMO

Enter:

    cd "$PDI/orbitals"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_homo.session.log

At the Multiwfn main menu, enter:

    5

Select:

    4

for orbital wavefunction.

When asked for the orbital index, enter:

    100

Select the medium-quality grid:

    2

When the post-processing menu appears, select:

    2

to export the grid as a Gaussian cube file.

Return to the main menu using the displayed return option, then enter:

    q

Inspect the generated file:

    ls -lt *.cub

Rename the newly generated orbital cube:

    mv <generated_HOMO_cube>.cub pdi_homo.cub

Replace `<generated_HOMO_cube>.cub` with the actual filename shown by `ls`.

### 4.2.2 Parent PDI LUMO

Run Multiwfn again:

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_lumo.session.log

Use:

    5
    4
    101
    2
    2

Return to the main menu and quit.

Rename the generated cube:

    mv <generated_LUMO_cube>.cub pdi_lumo.cub

### 4.2.3 Terminal-functionalized PDI HOMO

Enter:

    cd "$TFPDI/orbitals"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_homo.session.log

Use:

    5
    4
    140
    2
    2

Return to the main menu and quit.

Rename the generated cube:

    mv <generated_HOMO_cube>.cub \
      pdi_terminal_functionalized_homo.cub

### 4.2.4 Terminal-functionalized PDI LUMO

Run:

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_lumo.session.log

Use:

    5
    4
    141
    2
    2

Return to the main menu and quit.

Rename the generated cube:

    mv <generated_LUMO_cube>.cub \
      pdi_terminal_functionalized_lumo.cub

Create a record of the common visualization settings:

    cat > orbital_settings.txt <<'EOF'
    Orbitals: HOMO and LUMO
    Grid quality: Multiwfn medium-quality grid
    Recommended isovalue: ±0.03 a.u.
    Identical orientation, isovalue and rendering conventions must be used for both molecules.
    EOF

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

    cd "$PDI/charges"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_hirshfeld.session.log

At the main menu, enter:

    7

Select:

    1

for Hirshfeld atomic charges.

When asked how to obtain free-atom densities, select the option labelled:

    Use built-in sphericalized free-atom densities

When offered an export option, export the charges to a `.chg` file.

Return to the main menu and quit.

Inspect the generated files:

    ls -lt

Rename the charge file:

    mv <generated_charge_file>.chg pdi_hirshfeld.chg

### 4.3.2 Terminal-functionalized PDI

Enter:

    cd "$TFPDI/charges"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_hirshfeld.session.log

Use the same settings:

    7
    1

Export the charges to a `.chg` file.

Rename it:

    mv <generated_charge_file>.chg \
      pdi_terminal_functionalized_hirshfeld.chg

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

Avoid running games or other CPU-intensive applications during this calculation.

### 4.4.0 Test cube export before a long calculation

Before running a medium-quality ESP calculation, test whether the installed Multiwfn build can successfully export an ESP cube using a low-quality grid.

Enter:

    cd "$PDI/esp"

    export OMP_NUM_THREADS=8

    "$MWFN" "$PDI_WFN"

Use:

    5
    12
    2

This selects:

- grid-data analysis;
- total electrostatic potential;
- low-quality grid.

After the ESP calculation finishes, select:

    2

to export the grid as a Gaussian cube file.

Check whether a file such as the following was produced:

    totesp.cub

Exit using the displayed return option and then `q` from the main menu.

Check:

    ls -lh *.cub

If the low-quality export succeeds, delete the test file:

    rm -f totesp.cub

and proceed to the medium-quality calculation.

If selecting export option `2` produces:

    Fortran runtime error: Bad integer for item 1 in list input

do not immediately repeat the multi-hour calculation with the same executable. Test the official Multiwfn macOS binary using the low-quality grid first. Only run the medium-quality ESP calculation after cube export has been validated end-to-end.

### 4.4.1 Parent PDI electron-density cube

Enter:

    cd "$PDI/esp"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_density_cube.session.log

Use:

    5
    1
    2

This selects:

- grid-data analysis;
- electron density;
- medium-quality grid.

At the post-processing menu, select:

    2

to export a Gaussian cube file.

Return to the main menu and quit.

Rename the output:

    mv density.cub pdi_density.cub

Confirm:

    ls -lh pdi_density.cub

### 4.4.2 Parent PDI ESP cube

Enter:

    caffeinate -i env OMP_NUM_THREADS=8 \
      "$MWFN" "$PDI_WFN"

Use:

    5
    12
    2

This selects:

- grid-data analysis;
- total electrostatic potential;
- medium-quality grid.

Let the calculation finish completely.

At the post-processing menu, select:

    2

to export the ESP grid as a Gaussian cube file.

Return to the main menu and quit.

Rename the output:

    mv totesp.cub pdi_esp.cub

Confirm:

    ls -lh pdi_density.cub pdi_esp.cub

The density and ESP cubes should use the same grid quality and spatial region if they will be combined or compared point-by-point.

### 4.4.3 Parent PDI quantitative surface ESP

Run:

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_esp_surface.session.log

At the main menu, enter:

    12

Then enter:

    0
    
Start the surface analysis and wait for Multiwfn to finish evaluating the molecular surface.

The output should report quantities such as:

- global surface ESP minimum;

- global surface ESP maximum;

- local ESP extrema;

- molecular polarity index;

- positive and negative surface areas;

- polar and nonpolar surface fractions;

- average ESP and related surface statistics.

When the surface-analysis post-processing menu appears, enter:

    1

to export the extrema and surface statistics to:

    surfanalysis.txt

Then enter:

    2

to export the coordinates of the surface extrema to:

    surfanalysis.pdb

Return to the quantitative surface-analysis menu using:

    -1

Then return to the Multiwfn main menu using:

    -1

At the main menu, exit gracefully with:

    q

Back in Terminal, rename the exported files:

    mv surfanalysis.txt \

      pdi_esp_surface_statistics.txt

    mv surfanalysis.pdb \

      pdi_esp_extrema.pdb

Confirm that both files exist:

    ls -lh \

      pdi_esp_surface_statistics.txt \

      pdi_esp_extrema.pdb

The text file contains the numerical surface-ESP results, while the PDB file contains the positions of the ESP extrema for visualization in programs such as VMD, VESTA, ChimeraX, or Avogadro.

### 4.4.4 Terminal-functionalized PDI density cube

Enter:

    cd "$TFPDI/esp"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_density_cube.session.log

Use:

    5
    1
    2
    2

Return to the main menu and quit.

Rename:

    mv density.cub \
      pdi_terminal_functionalized_density.cub

### 4.4.5 Terminal-functionalized PDI ESP cube

Enter:

    caffeinate -i env OMP_NUM_THREADS=8 \
      "$MWFN" "$TFPDI_WFN"

Use:

    5
    12
    2

After the calculation finishes, export using:

    2

Return to the main menu and quit.

Rename:

    mv totesp.cub \
      pdi_terminal_functionalized_esp.cub

### 4.4.6 Terminal-functionalized PDI quantitative surface ESP

Run:

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_esp_surface.session.log

Use the same surface settings as for parent PDI:

    Electron-density isovalue: 0.001 a.u.
    Mapped function: Total electrostatic potential
    Grid spacing: approximately 0.20–0.25 Å

Start the surface analysis and wait for Multiwfn to finish evaluating the molecular surface.

The output should report quantities such as:

- global surface ESP minimum;

- global surface ESP maximum;

- local ESP extrema;

- molecular polarity index;

- positive and negative surface areas;

- polar and nonpolar surface fractions;

- average ESP and related surface statistics.

When the surface-analysis post-processing menu appears, enter:

    1

to export the extrema and surface statistics to:

    surfanalysis.txt

Then enter:

    2

to export the coordinates of the surface extrema to:

    surfanalysis.pdb

Return to the quantitative surface-analysis menu using:

    -1

Then return to the Multiwfn main menu using:

    -1

At the main menu, exit gracefully with:

    q

Back in Terminal, rename the exported files:

    mv surfanalysis.txt pdi_terminal_functionalizd_esp_surface_statistics.txt

    mv surfanalysis.pdb pdi_terminal_functionalizd_esp_extrema.pdb

Confirm that both files exist:

    ls -lh pdi_terminal_functionalizd_esp_surface_statistics.txt pdi_terminal_functionalizd_esp_extrema.pdb

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

    cd "$PDI/bond_orders"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_mayer.session.log

At the main menu, enter:

    9

Select:

    1

for Mayer bond-order analysis.

Export the complete Mayer bond-order matrix when offered.

Return to the main menu and quit.

Inspect:

    ls -lt

Rename the matrix:

    mv <generated_Mayer_file>.txt \
      pdi_mayer_bond_orders.txt

### 4.5.2 Terminal-functionalized PDI

Enter:

    cd "$TFPDI/bond_orders"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_mayer.session.log

Use:

    9
    1

Export the matrix.

Rename it:

    mv <generated_Mayer_file>.txt \
      pdi_terminal_functionalized_mayer_bond_orders.txt

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

    cd "$PDI/qtaim"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_qtaim.session.log

At the main menu, enter:

    2

Run the options labelled approximately:

    Search critical points starting from nuclear positions
    Search critical points starting from atom-pair midpoints
    Search critical points starting from triangle centres

For Multiwfn 3.8, these are commonly accessed using:

    2
    3
    4

Generate paths connecting nuclei and bond critical points using the option labelled:

    Generate paths connecting nuclei and bond critical points

This is commonly option:

    8

Use the critical-point management menu to export:

    CPs.txt
    CPs.pdb

Use the path management menu to export:

    paths.txt
    paths.pdb

Export properties evaluated at critical points:

    CPprop.txt

Rename:

    mv CPs.txt pdi_CPs.txt
    mv CPs.pdb pdi_CPs.pdb
    mv paths.txt pdi_paths.txt
    mv paths.pdb pdi_paths.pdb
    mv CPprop.txt pdi_CPprop.txt

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

    cd "$TFPDI/qtaim"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_qtaim.session.log

Repeat the same critical-point searches:

    2
    3
    4

Generate bond paths and export:

    CPs.txt
    CPs.pdb
    paths.txt
    paths.pdb
    CPprop.txt

Rename:

    mv CPs.txt \
      pdi_terminal_functionalized_CPs.txt

    mv CPs.pdb \
      pdi_terminal_functionalized_CPs.pdb

    mv paths.txt \
      pdi_terminal_functionalized_paths.txt

    mv paths.pdb \
      pdi_terminal_functionalized_paths.pdb

    mv CPprop.txt \
      pdi_terminal_functionalized_CPprop.txt

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

    cd "$PDI/elf_lol"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_elf.session.log

Use:

    5
    9
    2

At the post-processing menu, select:

    2

to export the cube file.

Return to the main menu and quit.

Inspect:

    ls -lt *.cub

Rename:

    mv <generated_ELF_cube>.cub pdi_elf.cub

### 4.7.2 Parent PDI LOL

Run:

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_lol.session.log

Use:

    5
    10
    2
    2

Return to the main menu and quit.

Rename:

    mv <generated_LOL_cube>.cub pdi_lol.cub

### 4.7.3 Terminal-functionalized PDI ELF

Enter:

    cd "$TFPDI/elf_lol"

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_elf.session.log

Use:

    5
    9
    2
    2

Rename:

    mv <generated_ELF_cube>.cub \
      pdi_terminal_functionalized_elf.cub

### 4.7.4 Terminal-functionalized PDI LOL

Run:

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_lol.session.log

Use:

    5
    10
    2
    2

Rename:

    mv <generated_LOL_cube>.cub \
      pdi_terminal_functionalized_lol.cub

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

    cd "$REPO"

Create:

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

Make it executable:

    chmod +x scripts/list_molden_atoms.py

Generate the atom-group lists:

    python3 scripts/list_molden_atoms.py "$PDI_WFN" \
      | tee "$PDI/dos/atom_groups.txt"

    python3 scripts/list_molden_atoms.py "$TFPDI_WFN" \
      | tee "$TFPDI/dos/atom_groups.txt"

Inspect:

    cat "$PDI/dos/atom_groups.txt"
    cat "$TFPDI/dos/atom_groups.txt"

Use these generated atom lists rather than manually guessing atom indices.

### 4.8.2 Parent PDI DOS and element PDOS

Enter:

    cd "$PDI/dos"

    "$MWFN" "$PDI_WFN" 2>&1 | tee pdi_dos.session.log

At the main menu, enter:

    10

In the DOS menu, configure:

    Energy unit: eV
    Energy range: -12 to 4 eV
    Broadening function: Gaussian
    FWHM: 0.30 eV

Open the fragment-definition menu.

Define:

    Fragment 1: all carbon atoms
    Fragment 2: all nitrogen atoms
    Fragment 3: all oxygen atoms
    Fragment 4: all hydrogen atoms

Copy the corresponding atom-index lists from:

    atom_groups.txt

For the initial demonstration, retain the default Mulliken orbital-composition method.

Generate the TDOS and fragment PDOS curves.

Export the curve data to text files.

Rename the outputs descriptively:

    pdi_tdos.txt
    pdi_pdos_C.txt
    pdi_pdos_N.txt
    pdi_pdos_O.txt
    pdi_pdos_H.txt

Record the settings:

    cat > pdi_dos_settings.txt <<'EOF'
    Energy unit: eV
    Energy range: -12 to 4 eV
    Broadening: Gaussian
    FWHM: 0.30 eV
    Projection fragments: C, N, O and H
    Orbital-composition method: Mulliken
    EOF

### 4.8.3 Terminal-functionalized PDI DOS and element PDOS

Enter:

    cd "$TFPDI/dos"

    cat atom_groups.txt

Run:

    "$MWFN" "$TFPDI_WFN" 2>&1 \
      | tee pdi_terminal_functionalized_dos.session.log

At the main menu, enter:

    10

Use exactly the same DOS settings:

    Energy unit: eV
    Energy range: -12 to 4 eV
    Broadening function: Gaussian
    FWHM: 0.30 eV
    Projection fragments: C, N, O and H
    Orbital-composition method: Mulliken

Export and rename:

    pdi_terminal_functionalized_tdos.txt
    pdi_terminal_functionalized_pdos_C.txt
    pdi_terminal_functionalized_pdos_N.txt
    pdi_terminal_functionalized_pdos_O.txt
    pdi_terminal_functionalized_pdos_H.txt

Because the functionalized molecule contains more atoms and molecular orbitals, its absolute TDOS intensity will naturally be larger.

For final comparison figures, either:

- normalize the TDOS curves;
- compare relative fragment contributions;
- compare orbital positions rather than absolute peak heights.

---

## 📁 4.9 Expected output structure

After completing Sections 4.1–4.8, each analysis folder should contain files resembling:

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
        ├── *_tdos.txt
        ├── *_pdos_*.txt
        ├── *_dos_settings.txt
        └── *_dos.session.log

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
8. TDOS and fragment PDOS data have been exported;
9. all comparison plots use consistent axes, colour scales, orientations and isovalues;
10. all menu selections, software versions and analysis settings are preserved in session logs or settings files.