import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType

from ase import Atoms, units
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
from ase.io.vasp import read_vasp
from ase.md.langevin import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary


import os
import shutil
import glob
from pathlib import Path
from time import perf_counter
from copy import deepcopy



print(f"pfp_api_client: {pfp_api_client.__version__}")

#############################Calculator setting###################################
calc_mode="CRYSTAL" # including +U correction
model_version="v4.0.0"  # the latest model version
method_type=EstimatorMethodType.PFVM
estimator = Estimator(calc_mode=calc_mode, model_version=model_version, method_type=method_type)
calculator = ASECalculator(estimator)
##################################################################################



def run_md_simulation_nvt(atoms, temperature, num_md_steps, num_interval, friction_coeff, 
                         traj_filename, log_filename, time_step):
    """Run NVT MD simulation with Berendsen thermostat"""
    
    # Initialize velocities according to Maxwell-Boltzmann distribution
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature, force_temp=True)
    Stationary(atoms)  # Set zero total momentum to avoid drifting
    
    # Set up Berendsen NVT dynamics
    taut = 1.0  # fs
    dyn = NVTBerendsen(atoms, time_step*units.fs, temperature_K=temperature, 
                       taut=taut*units.fs, loginterval=num_interval, 
                       trajectory=traj_filename)

    # Print statements
    def print_dyn():
        imd = dyn.get_number_of_steps()
        etot = atoms.get_total_energy()
        temp_K = atoms.get_temperature()
        elapsed_time = perf_counter() - start_time
        print(f"  {imd: >3}   {etot:.3f}    {temp_K:.2f}    {elapsed_time:.3f}")

    # Attach printing and logging
    dyn.attach(print_dyn, interval=num_interval)
    dyn.attach(MDLogger(dyn, atoms, log_filename, header=True, peratom=True, mode="a"), 
               interval=num_interval)

    # Run the dynamics
    start_time = perf_counter()
    print(f"    imd     Etot(eV)    T(K)    elapsed_time(sec)")
    dyn.run(int(num_md_steps))
    simulation_time = perf_counter() - start_time
    print(f"Total simulation time: {simulation_time:.3f} seconds")

def main():
    # Simulation parameters
    size_num = 2
    oxi = 30
    rand_num = 3
    file_names = ['37', '56', '73']
    
    # Load trajectory
    traj = Trajectory(f"../../3_O2_activation_input_generation/c_240722_new_pt_num/out_files/"
                     f"{file_names[size_num]}CA_{oxi}p_5rand.traj")
    atoms = traj[rand_num]
    name = file_names[size_num]
    atoms.calc = calculator
    
    # MD parameters
    time_step = 1    # fsec
    temperature = 500 + 273  # K
    num_md_steps = int(1E5)
    num_interval = int(1E3)
    friction_coeff = 0.002
    
    # Output files
    traj_filename = f"./240722_new_pt_{file_names[size_num]}CA_{oxi}p_{rand_num}_test.traj"
    log_filename = f"./240722_new_pt_{file_names[size_num]}CA_{oxi}p_{rand_num}_test.log"
    
    # Run MD simulation
    run_md_simulation_nvt(atoms, temperature, num_md_steps, num_interval, 
                         friction_coeff, traj_filename, log_filename, time_step)

if __name__ == "__main__":
    main() 