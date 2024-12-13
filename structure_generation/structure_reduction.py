import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType

from ase import Atoms
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

import pandas as pd
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
from collections import Counter

class StructureGenerator:
    def __init__(self):
        # Calculator setting
        self.calc_mode = "CRYSTAL"  # including +U correction
        self.model_version = "v4.0.0"  # the latest model version
        self.method_type = EstimatorMethodType.PFVM
        self.estimator = Estimator(
            calc_mode=self.calc_mode, 
            model_version=self.model_version, 
            method_type=self.method_type
        )
        self.calculator = ASECalculator(self.estimator)



    def load_structures(self, trajectory_path_template):
        """Load structures from trajectory files"""
        cluster_numbers = [0, 1, 2, 3, 4, 5]
        structures = []
        
        for i in cluster_numbers:
            traj_path = trajectory_path_template.format(i)
            log_filename = Trajectory(traj_path)
            optimized_str = log_filename[-1]  # Get last frame
            structures.append(optimized_str)
            
        return structures

    def create_ovacancy(self, cluster, vacancy_in_percent, O_consider_indices, center):

        # Get Ce and O indices
        ce_indices = np.array([atom.index for atom in cluster if atom.symbol == 'Ce'])
        o_indices = np.array([atom.index for atom in cluster if atom.symbol == 'O'])
        
        # Calculate oxygen numbers
        o_current = len(o_indices) - 3072  # Current O count
        o_ideal = len(ce_indices) * 2  # Ideal O count for CeO2 stoichiometry
        o_vac_num = int(len(ce_indices) * 2 * vacancy_in_percent / 100)  # Target vacancy count
        
        print(f"Current O: {o_current}, Ideal O: {o_ideal}, Target vacancies: {o_vac_num}")
        
        # Calculate number of O atoms to remove
        if o_ideal - o_current > o_vac_num:
            o_del_num = 0
        else:
            o_del_num = o_vac_num - (o_ideal - o_current)
        print(f"O atoms to remove: {o_del_num}")
        
        # Calculate distances and sort
        distances = [np.linalg.norm(cluster[o_index].position - center) 
                    for o_index in O_consider_indices]
        sorted_distances_indices = np.argsort(distances)[::-1]
        
        # Select indices to remove
        indices_to_remove_o_to_consider = sorted_distances_indices[:o_del_num]
        indices_to_remove = [O_consider_indices[idx] for idx in indices_to_remove_o_to_consider]
        
        print(f"Structure size before removal: {len(cluster)}")
        del cluster[indices_to_remove]
        print(f"Structure size after removal: {len(cluster)}")
        
        return cluster

    def generate_vacancy_structures(self, structures, O_indices_to_consider, center_list, 
                                  vacancy_percent=10):

        structures_with_vacancies = []
        
        for i, structure in enumerate(structures):
            # Create a copy of the structure to avoid modifying the original
            structure_copy = structure.copy()
            
            # Create vacancies
            modified_structure = self.create_ovacancy(
                structure_copy,
                vacancy_percent,
                O_indices_to_consider[i],
                center_list[i]
            )
            
            structures_with_vacancies.append(modified_structure)
            
        return structures_with_vacancies

    def select_tuning_O(self, cluster, cutoff_distance=3.0):
        """Select oxygen atoms for tuning based on distance criteria"""
        ce_indices = np.array([atom.index for atom in cluster if atom.symbol == 'Ce'])
        ce_positions = np.array([atom.position for atom in cluster if atom.symbol == 'Ce'])
        o_indices = np.array([atom.index for atom in cluster if atom.symbol == 'O'])
        
        ce_center_x = np.mean(ce_positions[:, 0])
        ce_center_y = np.mean(ce_positions[:, 1])
        min_z = np.min(ce_positions[:, 2])
        max_z = np.max(ce_positions[:, 2])
        ce_center = np.array([ce_center_x, ce_center_y, min_z])
        radius = max_z - min_z
        
        tuning_O_indices = []
        for o_index in o_indices:
            if cluster[o_index].position[2] - ce_center[2] > -3:
                closest_ce_distance = min(np.linalg.norm(cluster[o_index].position - cluster[ce_index].position) 
                                       for ce_index in ce_indices)
                if closest_ce_distance <= cutoff_distance:
                    tuning_O_indices.append(o_index)
        
        print(f"Number of O atoms selected for tuning: {len(tuning_O_indices)}")
        return tuning_O_indices

    def axis_dependent_size(self, atoms, axis):
        """Calculate size along specified axis"""
        ce_indices = np.array([atom.index for atom in atoms if atom.symbol == 'Ce'])
        ce_positions = np.array([atom.position for atom in atoms if atom.symbol == 'Ce'])
        ce_center_x = np.mean(ce_positions[:, 0])
        ce_center_y = np.mean(ce_positions[:, 1])
        min_z = np.min(ce_positions[:, 2])

        if axis == 'x':
            max_x = np.max(ce_positions[:, 1])
            radius = max_x - ce_center_x
        elif axis == 'y':
            max_y = np.max(ce_positions[:, 1])
            radius = max_y - ce_center_y
        else:
            max_z = np.max(ce_positions[:, 2])
            radius = max_z - min_z
            
        center = [ce_center_x, ce_center_y, min_z]
        return center, round(radius, 3)

    def print_cluster_info(self, atoms):
        """Print cluster information and return center"""
        center, size_x = self.axis_dependent_size(atoms, 'x')
        center, size_y = self.axis_dependent_size(atoms, 'y')
        center, size_z = self.axis_dependent_size(atoms, 'z')

        print(f"N_atoms: {len(atoms)}")
        print(f"Size_width_x: {round(size_x*2/10, 3)}")
        print(f"Size_width_y: {round(size_y*2/10, 3)}")
        print(f"Size_height: {round(size_z*2/10, 3)}")
        
        symbols = atoms.get_chemical_symbols()
        element_counts = Counter(symbols)
        element_counts_line = ', '.join(f"{element}: {count}" for element, count in element_counts.items())
        print(f"Composition: {element_counts_line}")
        print("===============================")
        return center

def main():
    generator = StructureGenerator()
    
    print(f"pfp_api_client: {pfp_api_client.__version__}")
    
    # Load structures
    trajectory_template = "../a_initial_cut/output/md_str_opt_20ps_{}.traj"
    structures = generator.load_structures(trajectory_template)
    
    if structures:
        # Generate O indices and centers for each structure
        O_indices_to_consider = []
        center_list = []
        
        for structure in structures:
            O_indices_to_consider.append(generator.select_tuning_O(structure))
            center = generator.print_cluster_info(structure)
            center_list.append(center)
        
        # Create structures with vacancies
        structures_with_vacancies = generator.generate_vacancy_structures(
            structures,
            O_indices_to_consider,
            center_list
        )
        
        # Save results
        if structures_with_vacancies:
            print(structures_with_vacancies[0].pbc)
            write('./input/2_10per_test.traj', structures_with_vacancies)

if __name__ == "__main__":
    main() 