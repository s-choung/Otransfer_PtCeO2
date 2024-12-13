# Otransfer_PtCeO<sub>2</sub>

This repository contains scripts for the preparation, simulation, and analysis of \(O_2\) uptake on Pt/CeO<sub>2</sub>-Al<sub>2</sub>O<sub>3</sub> systems.  
The work employs neural network potential (NNP) models to study oxygen activation and lattice oxygen transfer dynamics across varying ceria domain sizes.

## Files
- **Model Generation Scripts**: Scripts for generating Al<sub>2</sub>O<sub>3</sub> slabs, CeO<sub>2</sub> nanoparticles, and Pt-anchored configurations.
- **Simulation Setup**: Code for setting up \(O_2\) uptake simulations, including random \(O_2\) environment generation.
- **MD Simulation Scripts**: Scripts for running MD simulations with ASE and Matlantis.
- **Post-Simulation Analysis**: Tools for analyzing Ce-O coordination and visualizing results.

## Large Files
Due to GitHub's file size limitations, the following files are stored on Zenodo:
- Model structures for all domain sizes and Ce<sup>3+</sup> ratios.
- MD trajectories of annealing and \(O_2\) uptake simulations.

The Zenodo repository can be accessed at: [Zenodo link]

## Requirements
- Python 3.x
- ASE
- Matlantis (or other equivalent machine learning potential)
- Packmol
