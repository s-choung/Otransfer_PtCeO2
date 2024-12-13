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

## Structure Data Requirements
Each code file requires specific structure data to function correctly:

1. **`structure_generation/ceo2_al2o3_support.py`**
   - **Input:** Al<sub>2</sub>O<sub>3</sub> POSCAR file for slab generation.
   - **Output:** Generated Al<sub>2</sub>O<sub>3</sub> slab and CeO<sub>2</sub> clusters.

2. **`simulation/O2_md_simulation.py`**
   - **Input:** Trajectory files containing CeO<sub>2</sub> clusters with Pt atoms.
   - **Output:** MD simulation trajectory and log files.

3. **`structure_generation/o2_envronment_generation.py`**
   - **Input:** Structure files for Pt/CeO<sub>2</sub> systems and O<sub>2</sub> gas configurations.
   - **Output:** Packmol input files for generating O<sub>2</sub> environments.

4. **`structure_generation/pt_anchor.py`**
   - **Input:** Optimized structures of CeO<sub>2</sub> clusters.
   - **Output:** Pt-anchored structures saved in PDB format.

5. **`structure_generation/structure_reduction.py`**
   - **Input:** Trajectory files of optimized structures for vacancy generation.
   - **Output:** Structures with oxygen vacancies saved in trajectory format.

6. **`structure_generation/structure_generation.py`**
   - **Input:** Initial cluster structures for generating reduced structures.
   - **Output:** Reduced structures with specified vacancies.

7. **`structure_generation/pt_anchor.py`**
   - **Input:** Initial cluster structures for Pt anchoring.
   - **Output:** Structures with Pt atoms anchored to CeO<sub>2</sub> clusters.

## Notes
- Ensure that all required input files are available in the specified paths before running the scripts.
- The output files will be generated in the same directory or specified output directories as indicated in the scripts.
