from ase.cluster import Decahedron, Icosahedron, Octahedron, wulff_construction
from ase import Atoms 
from ase.io import Trajectory
import pandas as pd
from ase.build import bulk
from ase.constraints import ExpCellFilter, StrainFilter
from ase.optimize import LBFGS
import numpy as np
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType
from ase.io import Trajectory, write
from ase.build import bulk, surface, molecule, add_adsorbate, fcc111
from ase.constraints import ExpCellFilter, StrainFilter, FixAtoms, FixedPlane, FixBondLength
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.visualize import view
from ase.build.rotate import minimize_rotation_and_translation
from ase.io.vasp import read_vasp, write_vasp
from ase.io import read
import pandas as pd
import os
import shutil
import glob
from pathlib import Path
from PIL import Image, ImageDraw
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from scipy.spatial.distance import cdist
import subprocess

def pbc_cell(atoms, rod):
    """Apply periodic boundary conditions to atoms"""
    temp = atoms.copy()
    temp.cell = rod.cell
    temp.pbc = True
    return temp

def fix_atoms(atoms, rod):
    """Fix atoms below certain z-coordinate"""
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


class PackmolGenerator:
    def __init__(self):
        pass
        
    def generate_inp_file(self, rand_num, output_dir, structure_file, o2_gas, name, oxidation, o2_num, 
                         fixed_values=(38., 39., 8., 0., 0., 0.), 
                         box_dimensions=(3., 3., 10., 73., 75., 57.)):
        """
        Generate Packmol input file
        """
        content = f"""
tolerance 2
filetype pdb
seed -1

output ../out_files/{name}CA_{oxidation}p_{rand_num}_test.pdb
structure {structure_file}
  number 1
  center
  fixed {' '.join(map(str, fixed_values))}
end structure

structure {o2_gas}
   number {o2_num}           
   inside box {' '.join(map(str, box_dimensions))}
end structure
"""
        file_name = f"./{output_dir}/{name}CA_{oxidation}p_{rand_num}.inp"
        with open(file_name, 'w') as file:
            file.write(content.strip())
        return file_name

    def run_packmol(self, input_file):
        """Run Packmol with given input file"""
        input_dir = os.path.dirname(input_file)
        command = f"cd {input_dir} && packmol < {os.path.basename(input_file)}" # dependent on the packmol 
        
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Packmol ran successfully for {input_file}")
        else:
            print(f"Packmol encountered an error for {input_file}: {result.stderr}")

    def generate_multiple_structures(self, o2_nums, names, output_dir="inp_files", o2_gas='../O2_gas.pdb'):
        """Generate multiple structure files using Packmol"""
        for oxidation in [20, 30, 40]:
            for size in [0, 1, 2]:
                str_dir = f'./output/3_pt_anchor/{oxidation}p_new_Pt_{size}.pdb'
                for i in range(5):
                    rand_num = i
                    file_name = self.generate_inp_file(rand_num, output_dir, str_dir, o2_gas, names[size], 
                                                     oxidation, o2_nums[size])
                    self.run_packmol(file_name)

def main():
    # Initialize classes
    packmol = PackmolGenerator()
    
    # Example usage
    o2_nums = [152, 137, 119]
    names = [37, 56, 73]
    
    # Generate structures
    packmol.generate_multiple_structures(o2_nums, names)
    
    # Process generated structures
    ref_str_dir = './output/3_pt_anchor/40p_new_Pt_0.pdb'
    ref_str = read(ref_str_dir)
    
    for oxidation in [20, 30, 40]:
        for size in [0, 1, 2]:
            temp = []
            for i in range(5):
                str_dir = f'./output/4_o2_envronment_generation/{names[size]}CA_{oxidation}p_{i}.pdb'
                atoms_from_pdb_file = read(str_dir)
                atoms_from_pdb_file = fix_atoms(atoms_from_pdb_file, ref_str)
                atoms_from_pdb_file = pbc_cell(atoms_from_pdb_file, ref_str)
                temp.append(atoms_from_pdb_file)
            write(f'{names[size]}CA_{oxidation}p_5rand.traj', temp)

if __name__ == "__main__":
    main() 