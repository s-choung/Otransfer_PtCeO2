# Otransfer_PtCeO<sub>2</sub>

This repository contains scripts for the preparation, simulation, and analysis of O<sub>2</sub> uptake on Pt/CeO<sub>2</sub>-Al<sub>2</sub>O<sub>3</sub> systems.  
The work employs neural network potential (NNP) models to study oxygen activation and lattice oxygen transfer dynamics across varying ceria domain sizes.

## Repository Structure

### Structure Generation (`structure_generation/`)
- **Build Structures** (`build_structures.py`)
  - Al<sub>2</sub>O<sub>3</sub> slab generation
  - CeO<sub>2</sub> nanoparticle creation and hemisphere cleaving
  - Pt single atom anchoring on ceria domains
- **Geometry Optimization** (`optimize_geometry.py`)
  - O<sub>2</sub> adsorption optimization
  - Structure relaxation routines

### Simulation (`simulation/`)
- **MD Setup** (`md_setup.py`)
  - O<sub>2</sub> uptake simulation configuration
  - Random O<sub>2</sub> environment generation with Packmol
- **MD Execution** (`md_run.py`)
  - NVT ensemble implementation
  - Neural network potential integration via Matlantis

### Analysis (`analysis/`)
- **Structure Analysis** (`coordination.py`)
  - Ce-O coordination analysis
  - Surface/bulk Ce classification
- **Visualization** (`visualization.py`)
  - Lattice oxygen transfer visualization
  - O<sub>2</sub> activation pathway analysis

## Large Files
Due to GitHub's file size limitations, the following files are stored on Zenodo:
- Model structures for all domain sizes and Ce<sup>3+</sup> ratios
- MD trajectories of annealing and O<sub>2</sub> uptake simulations

The Zenodo repository can be accessed at: [Zenodo link]

## Requirements
- Python 3.x
- ASE (Atomic Simulation Environment)
- Matlantis (or other equivalent machine learning potential)
- Packmol
- NumPy
- Matplotlib
