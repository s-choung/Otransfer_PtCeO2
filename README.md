# Otransfer_PtCeO<sub>2</sub>

This repository contains scripts for the preparation, simulation, and analysis of O<sub>2</sub> uptake on Pt/CeO<sub>2</sub>-Al<sub>2</sub>O<sub>3</sub> systems.  
The work employs neural network potential (NNP) models to study oxygen activation and lattice oxygen transfer dynamics across varying ceria domain sizes.

### Structure Generation (`structure_generation/`)
- **Support Generation** (`ceo2_al2o3_support.py`)
  - Al<sub>2</sub>O<sub>3</sub> slab generation with vacuum spacing
  - CeO<sub>2</sub> cluster creation and hemisphere cutting
  - MD simulation for structure relaxation using NNP
- **Structure Reduction** (`structure_reduction.py`)
  - Oxygen vacancy creation based on distance criteria
  - Ce<sup>3+</sup>/Ce<sup>4+</sup> ratio control
- **Pt Anchoring** (`pt_anchor.py`)
  - Random Pt atom placement on CeO<sub>2</sub> domains
  - Energy-based structure optimization
- **O<sub>2</sub> Environment** (`o2_environment_generation.py`)
  - Random O<sub>2</sub> configuration generation using Packmol
  - System size-dependent O<sub>2</sub> number scaling

### Simulation (`simulation/`)
- **MD Simulation** (`O2_md_simulation.py`)
  - MD simulation of Pt/CeO<sub>2</sub>-Al<sub>2</sub>O<sub>3</sub> systems in O2 environment using NNP
  - NVT ensemble with Berendsen thermostat

### Analysis (`analysis/`)
- **Neighbor Analysis** (`1_neighbor_analysis.py`)
  - Oxygen neighbor identification
  - CeO<sub>2</sub> domain classification
- **Oxygen Activation Analysis** (`2_analyze_and_visualize.ipynb`)
  - Oxygen activation and lattice oxygen transfer analysis and visualization


## Requirements
- Python 3.x
- ASE (Atomic Simulation Environment)
- pfp_api_client (Neural Network Potential)
- Packmol
- NumPy
- Pymatgen
- SciPy

## Large Files

The following files are available through [Zenodo](https://zenodo.org/record/14498559):

- MD trajectories from annealing simulations (Geometry_relax_MD_CeO2_Al2O3.zip)
- O<sub>2</sub> uptake simulation trajectories of CeO2_Al2O3 (O2_MD_CeO2_Al2O3.zip)
- O<sub>2</sub> uptake simulation trajectories of Pt_CeO2_Al2O3 (O2_MD_Pt_CeO2_Al2O3.zip)


## Citation

Please cite the following paper when using this repository:

```
coming soon
```
