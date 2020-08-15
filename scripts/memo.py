from ase.io.vasp import read_vasp
from ase.io.xsd import write_xsd

contcar = read_vasp('CONTCAR')
write_xsd('contcar.xsd', contcar)

from pymatgen import Structure
struc = Structure.from_file('CONTCAR')

from pymatgen.io.vasp.inputs import Poscar

poscar = Poscar.from_file('CONTCAR')
poscar.structure





del atoms[atoms.numbers == 6]
del atoms[[atom.index for atom in atoms if atom.symbol == 'C']]
del atoms[[atom.symbol == 'C' for atom in atoms]]

del surf_2x2[[atom.position[2]/z_param < 0.35 for atom in surf_2x2]]
