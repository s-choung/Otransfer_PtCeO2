# Otransfer_PtCeO<sub>2</sub>

This repository contains scripts for the preparation, simulation, and analysis of O<sub>2</sub> uptake on Pt/CeO<sub>2</sub>-Al<sub>2</sub>O<sub>3</sub> systems.  
The work employs neural network potential (NNP) models to study oxygen activation and lattice oxygen transfer dynamics across varying ceria domain sizes.

## Repository Structure

### Structure Generation (`structure_generation/`)
- **Support Generation** (`ceo2_al2o3_support.py`)
  - Al<sub>2</sub>O<sub>3</sub> slab generation with vacuum spacing
  - CeO<sub>2</sub> cluster creation and hemisphere cutting
  - Structure optimization using NNP
  - MD simulation for structure relaxation
- **Pt Anchoring** (`pt_anchor.py`)
  - Random Pt atom placement on CeO<sub>2</sub> domains
  - Energy-based structure optimization
  - Multiple Pt configurations generation
- **Structure Reduction** (`structure_reduction.py`)
  - Oxygen vacancy creation based on distance criteria
  - Ce<sup>3+</sup>/Ce<sup>4+</sup> ratio control
  - Selective oxygen removal
- **O<sub>2</sub> Environment** (`o2_environment_generation.py`)
  - Packmol-based O<sub>2</sub> placement
  - Random O<sub>2</sub> configuration generation
  - System size-dependent O<sub>2</sub> number scaling

### Simulation (`simulation/`)
- **MD Simulation** (`O2_md_simulation.py`)
  - NVT ensemble with Berendsen thermostat
  - Neural network potential integration
  - Temperature and energy monitoring
  - Trajectory and log file generation

## Requirements
- Python 3.x
- ASE (Atomic Simulation Environment)
- pfp_api_client (Neural Network Potential)
- Packmol
- NumPy
- Pymatgen
- SciPy

## Structure Data Requirements

### Large Files
The following large files are required but not included in this repository:
- Model structures for all domain sizes and Ce<sup>3+</sup> ratios
- MD trajectories of annealing and O<sub>2</sub> uptake simulations

These files will be available through Zenodo after the publication of the work.

1. **`structure_generation/ceo2_al2o3_support.py`**
   - **Input:** Al<sub>2</sub>O<sub>3</sub> POSCAR file
   - **Output:** 
     - Al<sub>2</sub>O<sub>3</sub> slab with vacuum
     - Cut CeO<sub>2</sub> clusters
     - MD trajectories

2. **`structure_generation/pt_anchor.py`**
   - **Input:** Optimized CeO<sub>2</sub> structures
   - **Output:** 
     - Multiple Pt-anchored configurations
     - Energy-optimized structures
     - Trajectory files

3. **`structure_generation/structure_reduction.py`**
   - **Input:** Optimized structures
   - **Output:** Structures with controlled oxygen vacancies

4. **`structure_generation/o2_environment_generation.py`**
   - **Input:** 
     - Pt/CeO<sub>2</sub> structures
     - O<sub>2</sub> molecule template
   - **Output:** Packmol input files and configurations

5. **`simulation/O2_md_simulation.py`**
   - **Input:** Generated Pt/CeO<sub>2</sub>/O<sub>2</sub> structures
   - **Output:** 
     - MD trajectories
     - Energy/temperature logs

