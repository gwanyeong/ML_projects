# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:54:14 2020

@author: gyjung
"""

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun


with open('summary', 'w') as f:
    initial = Poscar.from_file('POSCAR', check_for_POTCAR = False)
    y_initial = initial.structure.lattice.abc[1]
    
    f.writelines('position\tlattice_y\tDisplacement_along_y\tfinal energy (eV)\n')
    for y_cnt in range(0,11):
        try:
            result = Vasprun('%02d_vertical/vasprun.xml' % (y_cnt))
            y_lattice = result.final_structure.lattice.abc[1]
            y_displacement = 0.5*(y_initial - y_lattice) # atom displacements = half of cell displacement
            final_E = result.final_energy
            if result.converged_electronic is False:
                f.writelines('%02d_vertical system is not electronically converged!' % (y_cnt))
            f.writelines('%02d\t%4.3f\t%4.3f\t%4.6f\n' % (y_cnt, y_lattice, y_displacement, final_E))
        except:
            f.writelines('%02d system is not properly calculated\n' % (y_cnt))
            continue

