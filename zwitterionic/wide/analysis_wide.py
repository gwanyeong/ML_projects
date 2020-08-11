# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:54:14 2020
@author: gyjung
"""

from pymatgen.io.vasp.outputs import Vasprun

with open('summary', 'w') as f:
    f.writelines('position\tDisplacement along a\tDisplacement along c\tfinal energy (eV)\n')
    for a_cnt in range(0,8):
        for c_cnt in range(0,8):
            try:
                result = Vasprun('%02d_%02d_shift/vasprun.xml' % (a_cnt, c_cnt))

                a_int = result.structures[-1].lattice.abc[0] * 0.1
                # b_int = result.structures[0].lattice.abc[1]
                c_int = result.structures[-1].lattice.abc[2] * 0.1

                final_E = result.final_energy
                if result.converged_electronic is False:
                    f.writelines('%02d_%02d system is not electronically converged!' % (a_cnt, c_cnt))
                f.writelines('(%02d,%02d)\t%4.3f\t%4.3f\t%4.6f\n' % (a_cnt, c_cnt, a_cnt * a_int, c_cnt * c_int, final_E))
            except:
                f.writelines('(%02d,%02d) system is not properly calculated\n' % (a_cnt, c_cnt))
                continue
~                                                                                                                                        
~                                                                                                                                        
~                                                                                                                                        
~                                                                                                                                        
~                                 
