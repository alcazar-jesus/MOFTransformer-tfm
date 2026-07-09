from ase.build import bulk
from ase.io import write, read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifParser
import numpy as np

# Create a test crystal
atoms = bulk('Cu', 'fcc', a=3.6)
write('test.cif', atoms)
write('test.extxyz', atoms)
write('test.xyz', atoms)

atoms_cif = read('test.cif')
atoms_extxyz = read('test.extxyz')
atoms_xyz = read('test.xyz')

print("CIF cell:", atoms_cif.cell.cellpar())
print("EXTXYZ cell:", atoms_extxyz.cell.cellpar())
print("XYZ cell:", atoms_xyz.cell.cellpar())
print("CIF pbc:", atoms_cif.pbc)
print("EXTXYZ pbc:", atoms_extxyz.pbc)
print("XYZ pbc:", atoms_xyz.pbc)
