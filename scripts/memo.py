from ase.io.vasp import read_vasp
from ase.io.xsd import write_xsd

contcar = read_vasp('CONTCAR')
write_xsd('contcar.xsd', contcar)

from pymatgen import Structure
struc = Structure.from_file('CONTCAR')

from pymatgen.io.vasp.inputs import Poscar

poscar = Poscar.from_file('CONTCAR')
poscar.structure
