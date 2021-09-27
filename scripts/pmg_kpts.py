from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath


struc = Structure.from_file('POSCAR')


x = HighSymmKpath(struc)

x.kpath
print(x.kpath['path'])
x.kpath['kpoints']

points = x.get_kpoints(line_density = 10) # default = 20

for i in range(len(points[0])):
    print(points[0][i],points[1][i])
