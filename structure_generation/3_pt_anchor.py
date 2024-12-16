import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType

from ase import Atoms, Atom
from ase.io import Trajectory, write
from ase.build import bulk, surface, molecule, add_adsorbate, fcc111
from ase.constraints import ExpCellFilter, StrainFilter, FixAtoms, FixedPlane, FixBondLength
from ase.optimize import LBFGS, BFGS, FIRE
from ase.neb import NEB
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.visualize import view
from ase.build.rotate import minimize_rotation_and_translation
from ase.md import MDLogger
from ase.io.vasp import read_vasp, write_vasp
from ase.io.proteindatabank import write_proteindatabank, read_proteindatabank

import pandas as pd
import ipywidgets as widgets
from IPython.display import display_png, Image as ImageWidget
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import os
import shutil
import glob
from pathlib import Path
from PIL import Image, ImageDraw

from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.vasp import Poscar
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

import numpy as np
from scipy.spatial.distance import cdist

class PtCeO2Generator:
    def __init__(self):
        print(f"pfp_api_client: {pfp_api_client.__version__}")
        
        # Calculator setting
        calc_mode = "CRYSTAL"  # including +U correction
        model_version = "v4.0.0"  # the latest model version
        method_type = EstimatorMethodType.PFVM
        self.estimator = Estimator(calc_mode=calc_mode, model_version=model_version, method_type=method_type)
        self.calculator = ASECalculator(self.estimator)

    def get_opt_energy(self, atoms, fmax=0.05, opt_mode: str = "normal"):    
        atoms.set_calculator(self.calculator)
        if opt_mode == "scale":
            opt1 = LBFGS(StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0]))
        elif opt_mode == "all":
            opt1 = LBFGS(ExpCellFilter(atoms))
        else:
            opt1 = LBFGS(atoms)
        opt1.run(fmax=fmax, steps=1000)
        return atoms.get_total_energy()

    def calculate_energy(self, atoms, fmax=0.05):
        atoms.calc = self.calculator
        E = self.get_opt_energy(atoms, fmax)
        return E


    def freeze(self, atoms):
        ce_indices = np.array([atom.index for atom in atoms if atom.symbol == 'Ce'])
        ce_positions = np.array([atom.position for atom in atoms if atom.symbol == 'Ce'])
        Al_positions = np.array([atom.position for atom in atoms if atom.symbol == 'Al'])
        Al_indices = np.array([atom.index for atom in atoms if atom.symbol == 'Al'])
        Al_z_axis_mean = np.mean(Al_positions[:, 2])
        freeze_indices = [atom.index for atom in atoms if atom.position[2] < Al_z_axis_mean+3]
        c = FixAtoms(indices=freeze_indices)
        print(len(freeze_indices))
        atoms.set_constraint(c)
        return atoms

    def Pt_anchor(self, cluster, pt_num):
        ce_indices = np.array([atom.index for atom in cluster if atom.symbol == 'Ce'])
        ce_positions = np.array([atom.position for atom in cluster if atom.symbol == 'Ce'])
        ce_center_x = np.mean(ce_positions[:, 0])
        ce_center_y = np.mean(ce_positions[:, 1])
        min_z, max_z = np.min(ce_positions[:, 2]), np.max(ce_positions[:, 2])
        center = [ce_center_x, ce_center_y, min_z]
        radius = max_z - min_z

        # Create a list to store the generated structures
        Pt_anchored_rand_structures = []

        for _ in range(10):
            # Create a copy of the original cluster
            new_cluster = cluster.copy()
            new_pt_positions = np.array([]).reshape(0, 3)

            # Generate and add multiple Pt atoms to the new cluster
            for _ in range(pt_num):
                while True:
                    r = radius + 2.5  # Distance from the center
                    theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
                    phi = np.random.uniform(0, np.pi)  # Polar angle
                    pt_x = center[0] + r * np.sin(phi) * np.cos(theta)
                    pt_y = center[1] + r * np.sin(phi) * np.sin(theta)
                    pt_z = center[2] + r * np.cos(phi)

                    # Check the distance to existing Pt atoms only when pt_num is greater than 1
                    if pt_num > 1 and np.all(cdist(new_pt_positions, [[pt_x, pt_y, pt_z]]) >= 8) and pt_z > min_z+3:
                        break
                    elif pt_num == 1 and pt_z > min_z+3:
                        break

                pt = Atom('Pt', position=[pt_x, pt_y, pt_z])
                new_cluster.append(pt)
                new_pt_positions = np.vstack([new_pt_positions, [pt_x, pt_y, pt_z]])

            # Append the new structure to the list
            Pt_anchored_rand_structures.append(new_cluster)

        return Pt_anchored_rand_structures

def main():
    generator = PtCeO2Generator()
    name = '40p' #Example with Ce3+ 40% (10% oxygen vacancy)
    
    # Load and prepare initial structures
    traj = Trajectory(f"./output/2_structure_reduction/CeO2_{name}.traj")
    opt_traj_list = []

    for atoms in traj:
        structure = atoms
        structure.set_pbc((True, True, True))
        structure.cell = [77.5509, 79.7665, 60]
        structure = generator.freeze(structure)
        print(structure.pbc)
        opt_traj_list.append(structure)

    # Fix periodic boundary conditions
    y_fixed_trajs_list = []
    for i, atoms in enumerate(opt_traj_list):
        x = atoms.cell[0][0]
        y = atoms.cell[1][1]
        z = atoms.cell[2][2]
        print(x, y, z)
        for atom in atoms:
            if atom.position[1] > y:
                atom.position = atom.position - [0, int(atom.position[1]/y)*y, 0]
            elif atom.position[1] < 0:       
                atom.position = atom.position - [0, (int(atom.position[1]/y)-1)*y, 0]
        y_fixed_trajs_list.append(atoms)
    print('Y_axis traj is done')

    x_y_fixed_trajs_list = []
    for i, atoms in enumerate(y_fixed_trajs_list):
        x = opt_traj_list[i].cell[0][0]
        y = opt_traj_list[i].cell[1][1]
        z = opt_traj_list[i].cell[2][2]
        for atom in atoms:
            if atom.position[0] > x:
                atom.position = atom.position - [int(atom.position[0]/x)*x, 0, 0]
            elif atom.position[0] < 0:       
                atom.position = atom.position - [(int(atom.position[0]/x)-1)*x, 0, 0]
        x_y_fixed_trajs_list.append(atoms)
    print(f'Xaxis_traj is done')
    print(len(x_y_fixed_trajs_list[1]))
    print(len(x_y_fixed_trajs_list[2]))

    # Generate Pt structures
    pt_num = [3, 9, 21]
    Pt_CeO2_supported_5per_rand_list = []
    for i, structure in enumerate(x_y_fixed_trajs_list):
        Pt_CeO2_supported_5per_0_rand = generator.Pt_anchor(structure, pt_num[i])
        Pt_CeO2_supported_5per_rand_list.append(Pt_CeO2_supported_5per_0_rand)

    # Save structures
    for i, structures in enumerate(Pt_CeO2_supported_5per_rand_list):
        write(f'./output/3_pt_anchor/{name}_new_Pt_rand_list_{i}.traj', structures)

    # Print center of mass for verification
    for i, structures in enumerate(Pt_CeO2_supported_5per_rand_list):
        for j, structure in enumerate(structures):
            print(structure.get_center_of_mass())

    # Calculate energies
    best_structure_list = []
    best_energy_list = []
    all_energy_list = []

    for i, structures in enumerate(Pt_CeO2_supported_5per_rand_list):
        energies = []
        for j, structure in enumerate(structures):
            e_tot = generator.calculate_energy(structure, fmax=0.05) 
            energies.append(e_tot)
        best_structure_index = np.argmin(energies)
        best_structure = Pt_CeO2_supported_5per_rand_list[i][best_structure_index]
        best_energy = energies[best_structure_index]
        best_structure_list.append(best_structure)
        best_energy_list.append(best_energy)
        all_energy_list.append(energies)
        write_proteindatabank(f'./output/3_pt_anchor/{name}_new_Pt_{i}.pdb',best_structure) 

if __name__ == "__main__":
    main() 