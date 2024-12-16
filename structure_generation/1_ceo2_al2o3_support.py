import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.io.vasp import read_vasp
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from collections import Counter
from copy import deepcopy
from ase.data import atomic_numbers, chemical_symbols, covalent_radii
from ase.md.langevin import Langevin
from ase.md import MDLogger
from ase import units, Atoms
from time import perf_counter

from pfp_api_client import __version__
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType
from ase.optimize import LBFGS
from ase.constraints import ExpCellFilter, StrainFilter

class CeO2Al2O3Generator:
    def __init__(self):
        # Initialize calculator
        calc_mode = "CRYSTAL"  # including +U correction
        model_version = "v4.0.0"  # the latest model version
        method_type = EstimatorMethodType.PFVM
        estimator = Estimator(calc_mode=calc_mode, model_version=model_version, method_type=method_type)
        self.calculator = ASECalculator(estimator)

    def get_opt_energy(self, atoms, fmax=0.001, opt_mode: str = "normal"):
        """Optimize structure and get energy"""
        atoms.set_calculator(self.calculator)
        if opt_mode == "scale":
            opt1 = LBFGS(StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0]))
        elif opt_mode == "all":
            opt1 = LBFGS(ExpCellFilter(atoms))
        else:
            opt1 = LBFGS(atoms)
        opt1.run(fmax=fmax)
        return atoms.get_total_energy()

    def generate_al2o3_slab(self, al2o3_poscar_path, vacuum_size=60.0):
        """Generate Al2O3 slab"""
        al2o3 = read_vasp(al2o3_poscar_path)
        
        slab_gen = SlabGenerator(
            initial_structure=AseAtomsAdaptor.get_structure(al2o3),
            miller_index=[0,0,1],
            min_slab_size=6.0,
            min_vacuum_size=vacuum_size,
            lll_reduce=False,
            center_slab=True,
            primitive=True,
            max_normal_search=1,
        )
        
        slabs = slab_gen.get_slabs(tol=0.3, bonds=None, max_broken_bonds=0, symmetrize=False)
        slab_atoms_list = [AseAtomsAdaptor.get_atoms(slab) for slab in slabs]
        slab = slab_atoms_list[0].copy()
        slab = slab * (16, 16, 1)
        
        # Shift slab to bottom of cell
        min_pos_z = np.min(slab.positions, axis=0)[2]
        slab.set_positions(slab.positions - [0, 0, min_pos_z])
        
        return slab

    def fix_slab_atoms(self, slab):
        """Fix atoms in the slab"""
        top = np.max(slab.positions[:, 2])
        mask = slab.positions[:, 2] < top + 0.1
        c = FixAtoms(mask=mask)
        slab.set_constraint(c)
        return slab, np.sum(mask)

    def get_xy_dimensions(self, atoms):
        """Get x,y dimensions of the slab"""
        return [round(atoms.cell[0, 0], 3), round(atoms.cell[1, 1], 3)]

    @staticmethod
    def size_based_on_center(atoms):
        """Calculate cluster size based on center of mass"""
        from scipy.spatial import distance
        center = atoms.get_center_of_mass()
        radius = distance.cdist([center], atoms.positions)
        diameter = round(np.max(radius)*2, 3)
        diameter_nm = diameter/10
        return round(diameter_nm, 3)

    @staticmethod
    def print_cluster_info(atoms, label):
        """Print information about atomic cluster"""
        size = CeO2Al2O3Generator.size_based_on_center(atoms)
        print(f"Size: {label}: {size} nm")
        print(f"N_atoms: {len(atoms)}")
        
        symbols = atoms.get_chemical_symbols()
        element_counts = Counter(symbols)
        element_counts_line = ', '.join(f"{element}: {count}" for element, count in element_counts.items())
        print(f"Composition: {element_counts_line}")
        print("===============================")

    def delete_half_depth(self, cluster):
        """Cut cluster in half based on position"""
        # Getting median position
        top = np.max(cluster.positions[:, 2])
        bottom = np.min(cluster.positions[:, 2])
        half_line = (top + bottom)/2
        
        # Delete the appropriate number of atoms from the bottom
        del_indices = [atom.index for atom in cluster if atom.position[2] < half_line]

        for index in reversed(del_indices):  # reverse to delete from the end
            del cluster[index]
        ce_indices = [atom.index for atom in cluster if atom.symbol == 'Ce']
        o_indices = [atom.index for atom in cluster if atom.symbol == 'O']
        indices = [ce_indices, o_indices]
        return cluster, indices

    def delete_half_atoms(self, cluster):
        """Cut cluster in half maintaining stoichiometry"""
        ce_indices = [atom.index for atom in cluster if atom.symbol == 'Ce']
        o_indices = [atom.index for atom in cluster if atom.symbol == 'O']
        
        # Sort these indices by z coordinate
        ce_indices.sort(key=lambda index: cluster[index].position[2])
        o_indices.sort(key=lambda index: cluster[index].position[2])

        # Calculate how many Ce and O atoms to delete to keep the 1:2 stoichiometry
        num_ce_to_delete = len(ce_indices) // 2
        num_o_to_delete = 2 * num_ce_to_delete

        # Collect indices to be deleted
        del_indices = ce_indices[:num_ce_to_delete] + o_indices[:num_o_to_delete]
        del_indices.sort()

        # Delete the appropriate number of atoms from the bottom
        for index in reversed(del_indices):
            del cluster[index]

        return cluster

    def add_cluster_to_support(self, slab, cluster):
        """Add cluster to Al2O3 support"""
        slab_sc = slab.copy()
        cluster = cluster.copy()

        slab_xy_size = np.min(slab_sc.cell.cellpar()[0:2])
        cluster_xy_size = np.max(
            (np.max(cluster.positions, axis=0) - np.min(cluster.positions, axis=0))[0:2]
        )

        slab_surface_xy_center = np.append(
            slab_sc.cell.cellpar()[0:2] / 2, np.max(slab_sc.positions, axis=0)[2]
        )

        cluster_surface_xy_center = np.append(
            np.mean(cluster.positions, axis=0)[0:2], np.min(cluster.positions, axis=0)[2]
        )

        cluster = Atoms(cluster.get_chemical_symbols(), 
                       cluster.positions - cluster_surface_xy_center)

        slab_surface_covalent_radii = covalent_radii[
            slab_sc.get_atomic_numbers()[np.argmax(slab_sc.positions, axis=0)[2]]
        ]

        cluster_surface_covalent_radii = covalent_radii[
            cluster.get_atomic_numbers()[np.argmin(cluster.positions, axis=0)[2]]
        ]

        cluster.translate(
            slab_surface_xy_center + [0, 0, slab_surface_covalent_radii + cluster_surface_covalent_radii]
        )

        supported = slab_sc.copy()
        supported += cluster
        
        return supported

    def run_md_simulation(self, atoms, temperature, num_md_steps, num_interval, 
                         friction_coeff, traj_filename, log_filename, time_step):
        """Run molecular dynamics simulation"""
        # Set up the MD simulation
        dyn = Langevin(atoms, time_step*units.fs, friction=friction_coeff, 
                       temperature_K=temperature, loginterval=num_interval, 
                       trajectory=traj_filename)

        # Print statements
        def print_dyn():
            imd = dyn.get_number_of_steps()
            etot = atoms.get_total_energy()
            temp_K = atoms.get_temperature()
            stress = atoms.get_stress(include_ideal_gas=True)/units.GPa
            stress_ave = (stress[0]+stress[1]+stress[2])/3.0 
            elapsed_time = perf_counter() - start_time
            print(f"  {imd: >3}   {etot:.3f}    {temp_K:.2f}    {stress_ave:.2f}  "
                  f"{stress[0]:.2f}  {stress[1]:.2f}  {stress[2]:.2f}  {stress[3]:.2f}  "
                  f"{stress[4]:.2f}  {stress[5]:.2f}    {elapsed_time:.3f}")

        dyn.attach(print_dyn, interval=num_interval)
        dyn.attach(MDLogger(dyn, atoms, log_filename, header=True, stress=True, 
                          peratom=False, mode="w"), interval=num_interval)

        # Run the dynamics
        start_time = perf_counter()
        print(f"    imd     Etot(eV)    T(K)    stress(mean,xx,yy,zz,yz,xz,xy)(GPa)  elapsed_time(sec)")
        dyn.run(num_md_steps)
        simulation_time = perf_counter() - start_time
        print(f"Total simulation time: {simulation_time:.3f} seconds")

def main():
    # Initialize generator
    generator = CeO2Al2O3Generator()
    
    # Generate Al2O3 slab
    al2o3_slab = generator.generate_al2o3_slab("./input/Al2O3.vasp")
    
    # Fix atoms
    fixed_slab, num_fixed = generator.fix_slab_atoms(al2o3_slab)
    print(f"Fixed atoms: {num_fixed} out of {fixed_slab.get_global_number_of_atoms()}")
    
    # Get dimensions
    dimensions = generator.get_xy_dimensions(fixed_slab)
    print(f"Slab dimensions: {dimensions}")
    
    # Example cluster sizes to analyze
    spherical_sizes = ['s','m','l']
    
    # Read clusters and create cut versions
    cut_spherical_clusters = []
    for size in spherical_sizes:
        cluster = read(f'./input/CeO2_{size}.xyz')
        cluster_copy = deepcopy(cluster)
        cut_cluster, indices = generator.delete_half_depth(cluster_copy)
        cut_spherical_clusters.append(cut_cluster)
        
        # Print cluster info
        print(f'Ce: {len(indices[0])}, O: {len(indices[1])}, '
              f'composition: {round(len(indices[0])/len(indices[1]),3)}')
        generator.print_cluster_info(cut_cluster, f'Cut {size}')

    # Add clusters to support
    supported_list = []
    for i, cluster in enumerate(cut_spherical_clusters):
        supported = generator.add_cluster_to_support(fixed_slab, cluster)
        supported_list.append(supported)
        write(f'./output/1_ceo2_al2o3/unrelaxed_ceo2_al2o3_{i}.vasp', supported)
    # Run MD simulations
    for i, structure in enumerate(supported_list):
        atoms = deepcopy(structure)
        atoms.pbc = True
        atoms.calc = generator.calculator
        
        # MD parameters
        time_step = 1.0    # fsec
        temperature = 773
        num_md_steps = 20000
        num_interval = 200
        friction_coeff = 0.002
        traj_filename = f"./output/md_str_opt_{i}.traj" # saved in in Zenodo with the name of Geometry_relax_MD_CeO2_Al2O3.zip
        log_filename = f"./output/md_relaxation_ceo2_al2o3_{i}.log"
        
        generator.run_md_simulation(atoms, temperature, num_md_steps, num_interval,
                                  friction_coeff, traj_filename, log_filename, time_step)
        write(f'./output/1_ceo2_al2o3/md_relaxed_ceo2_al2o3_{i}.vasp', atoms)  

if __name__ == "__main__":
    main() 