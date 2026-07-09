import os
import numpy as np
from ase.io import read
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor
from moftransformer.utils.prepare_data_xyz import get_crystal_graph

# Paths
dir_path = "/home/jalcazar/Documentos/Máster/TFM/prueba_cif_xyz/ejemplo_cif_xyz"
files = os.listdir(dir_path)
cif_file = [f for f in files if f.endswith('.cif')][0]
xyz_file = [f for f in files if f.endswith('.xyz') or f.endswith('.extxyz')][0]

cif_path = os.path.join(dir_path, cif_file)
xyz_path = os.path.join(dir_path, xyz_file)

# Read atoms
atoms_cif = read(cif_path)
atoms_xyz = read(xyz_path)

# Pymatgen structures
struct_cif = AseAtomsAdaptor.get_structure(atoms_cif)
struct_xyz = AseAtomsAdaptor.get_structure(atoms_xyz)

matcher = StructureMatcher()
is_match = matcher.fit(struct_cif, struct_xyz)

print(f"File CIF: {cif_file}")
print(f"File XYZ: {xyz_file}")
print(f"Structure Match: {is_match}")
print(f"Num Atoms: CIF={len(atoms_cif)} XYZ={len(atoms_xyz)}")
print(f"Cell CIF:\n{atoms_cif.cell[:]}")
print(f"Cell XYZ:\n{atoms_xyz.cell[:]}")

diff_pos = np.linalg.norm(atoms_cif.positions - atoms_xyz.positions)
print(f"L2 Norm of difference in positions (assuming same order): {diff_pos}")

# Compare graphs
atom_num_c, nbr_idx_c, nbr_dist_c, uni_idx_c, uni_count_c = get_crystal_graph(atoms_cif)
atom_num_x, nbr_idx_x, nbr_dist_x, uni_idx_x, uni_count_x = get_crystal_graph(atoms_xyz)

print("Atom nums equal:", np.array_equal(atom_num_c, atom_num_x))
if not np.array_equal(atom_num_c, atom_num_x):
    print("CIF atoms:", atom_num_c[:10])
    print("XYZ atoms:", atom_num_x[:10])

print("Nbr idx equal:", np.array_equal(nbr_idx_c, nbr_idx_x))
print("Nbr dist equal:", np.allclose(nbr_dist_c, nbr_dist_x, atol=1e-4))
print("Uni idx equal:", uni_idx_c == uni_idx_x)
