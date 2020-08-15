from ase.io.vasp import read_vasp
from ase.io.xsd import write_xsd

contcar = read_vasp('CONTCAR')
write_xsd('contcar.xsd', contcar)
