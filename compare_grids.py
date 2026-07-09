import os
import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cssr import Cssr
from moftransformer.utils.prepare_data_xyz import get_energy_grid
import pickle
import logging

dir_path = "/home/jalcazar/Documentos/Máster/TFM/prueba_cif_xyz/ejemplo_cif_xyz"
files = os.listdir(dir_path)
cif_file = [f for f in files if f.endswith('.cif')][0]
xyz_file = [f for f in files if f.endswith('.xyz') or f.endswith('.extxyz')][0]

cif_path = os.path.join(dir_path, cif_file)
xyz_path = os.path.join(dir_path, xyz_file)

atoms_cif = read(cif_path)
atoms_xyz = read(xyz_path)

logger = logging.getLogger("test")

# Generate grid for CIF
get_energy_grid(atoms_cif, "ejemplo_cif", dir_path, logger)
# Generate grid for XYZ
get_energy_grid(atoms_xyz, "ejemplo_xyz", dir_path, logger)

# Load grids
grid_cif = pickle.load(open(os.path.join(dir_path, "ejemplo_cif.griddata16"), "rb"))
grid_xyz = pickle.load(open(os.path.join(dir_path, "ejemplo_xyz.griddata16"), "rb"))

print("Grids equal:", np.array_equal(grid_cif, grid_xyz))
if not np.array_equal(grid_cif, grid_xyz):
    diff = np.abs(grid_cif.astype(np.float32) - grid_xyz.astype(np.float32))
    print("Max diff:", np.max(diff))
    print("Mean diff:", np.mean(diff))

# Compare CSSR writes
struct_cif = AseAtomsAdaptor.get_structure(atoms_cif)
struct_xyz = AseAtomsAdaptor.get_structure(atoms_xyz)

Cssr(struct_cif).write_file(os.path.join(dir_path, "cif_test.cssr"))
Cssr(struct_xyz).write_file(os.path.join(dir_path, "xyz_test.cssr"))

with open(os.path.join(dir_path, "cif_test.cssr")) as f1, open(os.path.join(dir_path, "xyz_test.cssr")) as f2:
    if f1.read() == f2.read():
        print("CSSR files are identical")
    else:
        print("CSSR files are DIFFERENT")

