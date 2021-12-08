from ase.calculators.vasp import Vasp
from ase.dft.bandgap import bandgap
from ase.io import read

model = read('POSCAR')
calc = Vasp(restart = True)
model.calc = calc
gap, p1, p2 = bandgap(model.calc)
print(gap," ",p1," ",p2)
