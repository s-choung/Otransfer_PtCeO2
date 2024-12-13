# Otransfer_PtCeO2

## Description
This repository contains scripts for the preparation, simulation, and analysis of O2 uptake on Pt/CeO2-Al2O3 systems. 
The work employs neural network potential (NNP) models to study oxygen activation and lattice oxygen transfer dynamics across varying ceria domain sizes.

## Files
- **Model Generation Scripts**: Scripts for generating Al2O3 slabs, CeO2 nanoparticles, and Pt-anchored configurations.
- **Simulation Setup**: Code for setting up O2 uptake simulations, including random O2 environment generation.
- **MD Simulation Scripts**: Scripts for running MD simulations with ASE and Matlantis.
- **Post-Simulation Analysis**: Tools for analyzing Ce-O coordination and visualizing results.
- **Energy Minimization**: Geometry optimization and NEB setup scripts.

## Large Files
Due to GitHub's file size limitations, the following files are stored on Zenodo:
- Model structures for all domain sizes and Ce3+ ratios.
- MD trajectories of annealing and O2 uptake simulations.
- Optimized structures of oxygen adsorption and NEB results.

The Zenodo repository can be accessed at: [Zenodo link]

## Requirements
- Python 3.x
- ASE
- Matlantis (or other equvalent machine learning potential)
- Packmol
