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