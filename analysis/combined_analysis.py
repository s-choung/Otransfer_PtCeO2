# Essential imports
import numpy as np
from ase.io import Trajectory
from ase.constraints import FixAtoms
from scipy.spatial import distance, cKDTree
from tqdm import tqdm
import os
import pickle  # This was missing but is needed for pickle.dump operations

# Global settings
CUTOFF_VAL = 3.5  # Ce-O cutoff distance


class TrajectoryAnalyzer:
    """Analyze atomic trajectories and identify important atomic indices"""
    def __init__(self, traj_list):
        self.traj_list = traj_list
        self.initial_structures = []
        self.initial_ce_centers = []
        self.initial_radiuses = []
        self.O_indices_to_consider = []
        self.initial_gas_o2_indices = []
        self.initial_Ce_indices = []
        self.initial_O_indices_to_consider_sorted_list = []
        self.initial_Ce_indices_sorted_list = []

    def calculate_ce_geometry(self, traj):
        ce_positions = np.array([atom.position for atom in traj if atom.symbol == 'Ce'])
        ce_indice = np.array([atom.index for atom in traj if atom.symbol == 'Ce'])
        ce_center_x, ce_center_y = np.mean(ce_positions[:, 0]), np.mean(ce_positions[:, 1])
        ce_min_z, ce_max_z = np.min(ce_positions[:, 2]), np.max(ce_positions[:, 2])
        ce_center = np.array([ce_center_x, ce_center_y, ce_min_z])
        distances = np.linalg.norm(ce_positions[:, np.newaxis] - ce_center, axis=2)
        radius = ce_max_z - ce_min_z
        return ce_center, radius, ce_max_z, ce_min_z, ce_indice

    def select_tuning_O(self, cluster, cutoff_distance=3.0):
        ce_center, radius, ce_max_z, ce_min_z, ce_indice = self.calculate_ce_geometry(cluster)
        o_indices_tot = np.array([atom.index for atom in cluster if atom.symbol == 'O' and atom.position[2] < ce_max_z + 3])
        o_indices_lower_than_ce_max = np.array([atom.index for atom in cluster if atom.symbol == 'O' and atom.position[2] < ce_max_z + 3])
        tuning_O_indices = []
        for o_index in o_indices_lower_than_ce_max:
            if cluster[o_index].position[2] - ce_center[2] > -2 and np.abs(cluster[o_index].position[0] - ce_center[0]) < radius and np.abs(cluster[o_index].position[1] - ce_center[1]) < radius:
                tuning_O_indices.append(o_index)
        print('number of O that is part of CeO2', len(tuning_O_indices))
        return tuning_O_indices

    def initial_gas_o2_indices_generator(self, cluster, cutoff_distance=3.0):
        ce_center, radius, ce_max_z, ce_min_z, ce_indice = self.calculate_ce_geometry(cluster)
        o_indices_tot = np.array([atom.index for atom in cluster if atom.symbol == 'O' and atom.position[2] > ce_min_z + 2])
        initial_gas_o2_indices = np.array([index for index in o_indices_tot if distance.euclidean(cluster[index].position, ce_center) > radius + 4])
        susspicious_o_indice = []
        for o2_gas_index in initial_gas_o2_indices:
            for ce_index in ce_indice:
                ce_position = cluster.positions[ce_index]
                if np.all(np.abs(cluster.positions[o2_gas_index] - ce_position) < 3):
                    susspicious_o_indice.append(o2_gas_index)
        print('number of O in O2 gas', len(initial_gas_o2_indices))
        return initial_gas_o2_indices

    def get_indices_by_distance_from_center(self, atoms, indices, center):
        index_distance_pairs = []
        for index in indices:
            dist = round(distance.euclidean(atoms[index].position, center), 3)
            index_distance_pairs.append((index, dist))
        sorted_pairs = sorted(index_distance_pairs, key=lambda x: x[1])
        sorted_indices = [pair[0] for pair in sorted_pairs]
        return sorted_indices

    def process_trajectories(self):
        for traj in self.traj_list:
            self.initial_structures.append(traj[0])
            ce_center, radius, ce_max_z, ce_min_z, ce_indice = self.calculate_ce_geometry(traj[0])
            self.initial_ce_centers.append(ce_center)
            self.initial_radiuses.append(radius)
            self.O_indices_to_consider.append(self.select_tuning_O(traj[0]))
            self.initial_gas_o2_indices.append(self.initial_gas_o2_indices_generator(traj[0]))
            self.initial_O_indices_to_consider_sorted_list.append(
                self.get_indices_by_distance_from_center(traj[0], self.O_indices_to_consider[-1], self.initial_ce_centers[-1]))
            self.initial_Ce_indices.append(np.array([atom.index for atom in traj[0] if atom.symbol == 'Ce']))
            self.initial_Ce_indices_sorted_list.append(
                self.get_indices_by_distance_from_center(traj[0], self.initial_Ce_indices[-1], self.initial_ce_centers[-1]))

class CeOBondsAnalyzer:
    """Analyze Ce-O bonds and their properties"""
    def __init__(self, traj_list, ce_centers):
        self.traj_list = traj_list
        self.ce_centers = ce_centers

    def get_indices_by_distance_from_center(self, atoms, center):
        ce_indices = np.where(atoms.symbols == 'Ce')[0]
        ce_positions = atoms.positions[ce_indices]
        distances = np.linalg.norm(ce_positions - center, axis=1)
        sorted_indices = ce_indices[np.argsort(distances)]
        return sorted_indices

    def near_Ce_information(self, Atoms, ce_indices):
        neighbor_o_indices_num_list = []
        neighbor_o_indices_list = []

        positions = Atoms.positions
        symbols = Atoms.symbols
        ce_positions = positions[ce_indices]
        o_positions = positions[symbols == 'O']
        tree = cKDTree(o_positions)
        Ce_dictionary = {}

        for i, (ce_index, ce_pos) in enumerate(zip(ce_indices, ce_positions)):
            neighbor_o_indices = tree.query_ball_point(ce_pos, CUTOFF_VAL)
            bond_num = len(neighbor_o_indices)
            neighbor_o_indices_num_list.append(bond_num)
            neighbor_o_indices_list.append(neighbor_o_indices)
            Ce_dictionary[i] = {
                "ce": np.round(ce_pos, 2),
                "ce_index": ce_index,
                "nei": bond_num,
                "nei_o_index": neighbor_o_indices
            }

        return Ce_dictionary

    @staticmethod
    def count_values_in_list(value_list, values_to_count):
        counts = {}
        for value in values_to_count:
            counts[value] = np.count_nonzero(value_list == value)
        return counts

    def process_trajectories(self):
        Ce_dictionary_traj_list = []
        for j, traj in enumerate(self.traj_list):
            Ce_dictionary_traj = []
            for atoms in tqdm(traj, desc="Processing"):
                ce_indices = self.get_indices_by_distance_from_center(atoms, self.ce_centers[j])
                Ce_dictionary = self.near_Ce_information(atoms, ce_indices)
                Ce_dictionary_traj.append(Ce_dictionary)
            Ce_dictionary_traj_list.append(Ce_dictionary_traj)
        return Ce_dictionary_traj_list

def apply_pbc_to_traj_list(traj_list, x, y):
    """Apply periodic boundary conditions to trajectory list"""
    pbc_traj_list = []
    for traj in traj_list:
        pbc_traj = []
        for i, atoms in enumerate(traj):
            for atom in atoms:
                # Apply PBC in x direction
                if atom.position[0] > x:
                    atom.position[0] -= x
                elif atom.position[0] < 0:
                    atom.position[0] += x
                    
                # Apply PBC in y direction
                if atom.position[1] > y:
                    atom.position[1] -= y
                elif atom.position[1] < 0:
                    atom.position[1] += y
            pbc_traj.append(atoms)
            if i == len(traj)-1:
                print(f"Y_axis traj {len(pbc_traj_list)} is done")
        pbc_traj_list.append(pbc_traj)
        print(f"X-axis traj {len(pbc_traj_list)-1} is done")
    return pbc_traj_list

def main():
    """Main execution function that combines functionality from both notebooks"""
    # Configuration
    list_item = [0, 1, 2]
    size = [37, 56, 73]
    titles = ['1_5per']#, '2_10per', '3_ce3p_30']
    rand_list = [0, 1, 2, 3, 4]

    # Process each title
    for title in titles:
        # Set oxidation state based on title
        if title == '1_5per':
            oxid = 20
        elif title == '2_10per':
            oxid = 40
        else:
            oxid = 30

        # Create trajectory lists
        traj_list_list = []
        for j in rand_list:
            traj_list = []
            for i in list_item:
                traj_filename = f"../../calculations/3_20240722_new_pt_num/done/240722_new_pt_{size[i]}CA_{oxid}p_{j}.traj"
                try:
                    traj = Trajectory(traj_filename)
                    print(f'{size[i]}CA_{oxid}p_{j}')
                    
                    if len(traj) < 100:
                        print(f'{size[i]}CA_{oxid}p_{j} traj short!!')
                        continue

                    every_5ps_traj = [traj[i] for i in list(range(0, 101, 5))]
                    traj_list.append(every_5ps_traj)
                except Exception as e:
                    print(f"Error processing trajectory {traj_filename}: {str(e)}")
                    continue
            
            if traj_list:  # Only append if traj_list is not empty
                traj_list_list.append(traj_list)

        if not traj_list_list:  # Skip if no trajectories were processed
            continue

        # Get cell parameters from first trajectory
        x = traj_list_list[0][0][0].cell[0][0]
        y = traj_list_list[0][0][0].cell[1][1]
        z = traj_list_list[0][0][0].cell[2][2]
        print(x, y, z)

        # Apply PBC and save results
        pbc_nested_traj = []
        for nth_traj_list in traj_list_list:
            pbc_traj_list = apply_pbc_to_traj_list(nth_traj_list, x, y)
            pbc_nested_traj.append(pbc_traj_list)

        # Save PBC trajectories
        with open(f'{title}_pbc_test.pkl', 'wb') as f:
            pickle.dump(pbc_nested_traj, f)

        # Process trajectories with analyzers
        analyzer_list = []
        for pbc_traj_list in pbc_nested_traj:
            analyzer = TrajectoryAnalyzer(pbc_traj_list)
            analyzer.process_trajectories()
            analyzer_list.append(analyzer)

        # Analyze Ce-O bonds
        Ce_dictionary_traj_list_list = []
        for i, pbc_traj_list in enumerate(pbc_nested_traj):
            print('<<<<<<<<<<', i, '>>>>>>>>>>>')
            ceo_analyzer_temp = CeOBondsAnalyzer(pbc_traj_list, analyzer_list[i].initial_ce_centers)
            Ce_dictionary_traj_list = ceo_analyzer_temp.process_trajectories()
            Ce_dictionary_traj_list_list.append(Ce_dictionary_traj_list)

        # Save results
        file_name = int(CUTOFF_VAL * 10)
        os.makedirs('./nei_test', exist_ok=True)
        with open(f'./nei_test/{title}_nei_{file_name}_test.pkl', 'wb') as f:
            pickle.dump(Ce_dictionary_traj_list_list, f)

if __name__ == "__main__":
    main()
