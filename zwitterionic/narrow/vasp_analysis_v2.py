# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:54:14 2020

@author: gyjung
"""

from pymatgen.io.vasp.outputs import Vasprun

with open('summary', 'w') as f:
    f.writelines('position\tDisplacement along x\tDisplacement along z\tfinal energy (eV)\n')
    for x_cnt in range(0,6):
        for z_cnt in range(0,6):
            try:
                result = Vasprun('%d_%d_shift/vasprun.xml' % (x_cnt, z_cnt))
                
                x_int = result.structures[-1].lattice.abc[0] * 0.01
                # y_int = result.structures[0].lattice.abc[1]
                z_int = result.structures[-1].lattice.abc[2] * 0.01

                final_E = result.final_energy
                if result.converged_electronic is False:
                    f.writelines('%d_%d system is not electronically converged!' % (x_cnt, z_cnt))
                f.writelines('(%d,%d)\t%4.3f\t%4.3f\t%4.6f\n' % (x_cnt, z_cnt, x_cnt * x_int, z_cnt * z_int, final_E))
            except:
                f.writelines('(%d,%d) system is not properly calculated\n' % (x_cnt, z_cnt))
                continue




