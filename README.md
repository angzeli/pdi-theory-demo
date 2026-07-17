# 🧪 PDI Theory Demo

A lightweight computational chemistry playground for developing and validating ORCA–Multiwfn–Python workflows on molecular photocatalysts.

This repository serves as a sandbox for learning quantum chemical calculations, automating post-processing, and producing publication-quality theoretical figures. The long-term goal is to establish a reproducible workflow that can be readily adapted to real research projects involving photocatalytic materials.

---

## 🎯 Objectives

This project aims to:

- Learn molecular DFT using ORCA
- Explore wavefunction analysis with Multiwfn
- Automate data extraction using Python
- Produce publication-quality figures
- Develop reproducible computational chemistry workflows
- Prepare for real-world photocatalysis research

---

## 🛠️ Planned Workflow

```text
    ChemDraw & Avogadro
        │
        3D Molecular Structure
        │
        ▼
    ORCA ground-state calculations
        │
        ├── Geometry optimization
        ├── Frequency validation
        └── Higher-level single-point calculation
        │
        ▼
    Ground-state Multiwfn
        │
        ├── HOMO / LUMO
        ├── Hirshfeld charges       
        ├── ESP
        ├── Mayer bond orders
        ├── QTAIM topology
        ├── ELF and LOL
        └── DOS, PDOS, OPDOS
        │
        ▼
    Python, VMD, IQmol, Powerpoint post-processing
        │
        ├── Data parsing
        ├── Plotting
        ├── Figure generation and combination
        └── Report automation
        │
        ▼
    ORCA follow-up calculations
        │
        ├── TDA
        ├── TD-DFT Analysis
        ├── Charged states
        ├── Adsorption complexes
        ├── Reaction intermediates
        └── Transition states where justified
        │
        ▼
    Multiwfn follow-up analysis
        │
        ├── NTO analysis
        ├── Hole–electron analysis
        ├── Fukui functions
        ├── NCI / IGMH
        ├── Charge-density difference
        └── Fragment charge transfer
        │
        ▼
    Python, VMD, IQmol, Powerpoint post-processing
        │
        ├── Data parsing
        ├── Plotting
        ├── Figure generation and combination
        └── Report automation
```

---

## 📂 Repository Structure

```text
pdi-theory-demo/
├── calculations/
├── structures/
├── scripts/
├── analysis/
├── figures/
├── results/
├── docs/
├── README.md
└── LICENSE
```

## 🔬 Planned Demonstration

The initial demonstration will compare:

- Parent PDI
- Terminally modified PDI

Possible analyses include:

- Geometry optimization
- Frontier molecular orbitals
- HOMO–LUMO gap
- Electrostatic potential (ESP)
- Atomic charge analysis
- O2 adsorption
- Charge transfer
- TD-DFT excitation
- Publication-ready visualizations

---

## 🚀 Tech Stack

- ORCA
- Multiwfn
- Python
- NumPy
- Matplotlib
- Jupyter Notebook (optional)
- ASE (planned)

---

## 📈 Future Plans

Potential future additions include:

- Automated ORCA input generation
- Automatic convergence checking
- Batch calculation management
- TD-DFT spectrum generation
- Publication-ready plotting utilities
- Composite figure generation
- Computational chemistry workflow automation

---

## 📚 Purpose

This repository is primarily intended as a learning and development environment for computational chemistry.

The code is designed to be modular, reproducible, and extensible so that components developed here can later be incorporated into larger research workflows and future scientific software projects.

---

## 👤 Author

Angze Li

MSci Chemistry
Imperial College London