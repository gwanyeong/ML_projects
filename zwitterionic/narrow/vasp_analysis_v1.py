# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:54:14 2020

@author: gyjung
"""

from pymatgen.io.vasp.outputs import Vasprun

with open('summary', 'w') as f:
    for x_int in range(0,6):
        for z_int in range(0,6):
            result = Vasprun('%d_%d_shift/vasprun.xml' % (x_int, z_int))
            final_E = result.final_energy
            f.writelines('(%d,%d)\t%4.6f\n' % (x_int, z_int, final_E))
            
