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