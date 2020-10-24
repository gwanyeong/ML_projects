# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 11:40:47 2020
@author: gyjung
"""

import os
import shutil
import itertools
import fileinput
import time
import warnings
warnings.filterwarnings('ignore')

from ase.io.vasp import read_vasp, write_vasp
from ase.build import surface, make_supercell
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.build import sort        

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

##############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

###############################################################################
LDH_list = ["ScH2O2", "TiH2O2","VH2O2", "CrH2O2","MnH2O2","FeH2O2","CoH2O2","NiH2O2","CuH2O2","ZnH2O2",
            "YH2O2","ZrH2O2","NbH2O2","MoH2O2","TcH2O2","RuH2O2","RhH2O2","PdH2O2","AgH2O2","CdH2O2",
            "HfH2O2","TaH2O2","WH2O2","ReH2O2","OsH2O2","IrH2O2","PtH2O2","AuH2O2","HgH2O2"]


###############################################################################
createFolder('models_slab')

###############################################################################
with open('models_slab/slab_modeling.log', 'w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_slab')  

    for idx, formula in enumerate(LDH_list):
        f.writelines('%02d_%s\n' % (idx + 1.0, formula))

        TM = formula.replace("H2O2","")
        file_path = '%02d_%s/opt/' % (idx + 1.0, formula)
        createFolder(file_path + 'slab_100')

        """
        Slab modeling
        """    
        bulk_opt = read_vasp(file_path + 'CONTCAR')
       
        slab = surface(bulk_opt, (1,0,0), 5, vacuum = 7.5)  
        slab_s = make_supercell(slab, [[2,0,0],[0,1,0],[0,0,1]]) 

        f.writelines('\tN_atoms(before): %d' % len(slab_s))
        
        #Slab modification
        z_list = []
        for z in range(len(slab_s)):
            z_list.append(slab_s.get_positions()[z,2])
        z_list.sort()

        z_lower = z_list[5]
        z_upper = z_list[-1]
        
        del slab_s[[atom.position[2] < z_lower + 0.01 for atom in slab_s]]
        del slab_s[[atom.position[2] > z_upper - 0.01 for atom in slab_s]]

        f.writelines('    N_atoms(after) : %d' % (len(slab_s)))
        
        #Layer fix
        fix_atom_list = []
        for atom in slab_s:
            if atom.position[2] / slab_s.cell[2,2] < 0.5:
                fix_atom_list.append(atom.index)

        slab_s.set_constraint(FixAtoms(indices = fix_atom_list))
        
        f.writelines('    N_atoms(fixed) : %d\n' % (len(fix_atom_list)))
 
#        view(slab_s)      
        
        slab_sorted = sort(slab_s)
        
        """
        INCAR
        """
        INCAR = Incar.from_file(file_path + 'INCAR')
        
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = 0
        INCAR['IDIPOL'] = 3
        INCAR['ICHARG'] = 2
        INCAR['ISPIN'] = 1
        INCAR['NSW'] = 0
        INCAR['IBRION'] = -1
        INCAR['LREAL'] = '.FALSE.'
        
        # Hubbard U setting
        elements = []
        for element in bulk_opt.get_chemical_symbols():
            if element not in elements:
                elements.append(element)
        
        # Correction for previous error
        if INCAR['LDAUL'] == [0.0, 0.0, 0.0]:
            INCAR['LDAUL'] = [-1, -1, -1]
    
        
        LDAUL_dic = dict(zip(elements, INCAR['LDAUL']))
        LDAUU_dic = dict(zip(elements, INCAR['LDAUU']))
        LDAUJ_dic = dict(zip(elements, INCAR['LDAUJ']))
       
        INCAR['LDAUL'] = [LDAUL_dic[x] for x in sorted(elements)]
        INCAR['LDAUU'] = [LDAUU_dic[x] for x in sorted(elements)]
        INCAR['LDAUJ'] = [LDAUJ_dic[x] for x in sorted(elements)]           
        
        f.writelines(['\t','elements: ',str(sorted(elements)),'\n'])

        """
        MAGMOM
        """
        mag_list = []
        for e in elements:
            if e == TM:
                mag_list.append(5.0)
            else:
                mag_list.append(0.6)

        mag_dic = dict(zip(elements, mag_list))
        
        magmom = []
        
        for el in slab_sorted.get_chemical_symbols():
            magmom.append(mag_dic[el])    
        
        INCAR['MAGMOM'] = magmom

        magm = []
        for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
            magm.append("{}*{}".format(len(tuple(g)), m))
            
        
        """
        POTCAR
        """

        potcar_ori = Potcar.from_file(file_path + 'POTCAR')
        potcar_symbols = potcar_ori.as_dict()['symbols']
        potcar_symbols.sort()
        POTCAR = Potcar(potcar_symbols) 
        f.writelines(['\t','POTCAR_symbols: ', str(potcar_symbols), '\n'])

        f.writelines(['\t', 'LDAUL: ',str(INCAR['LDAUL']),'\t',
                      'LDAUU: ',str(INCAR['LDAUU']),'\t',
                      'LDAUJ: ',str(INCAR['LDAUJ']),'\n'])

        f.writelines(['\t','MAGMOM: ',str(magm),'\n'])

        print(sorted(elements)," ",INCAR['LDAUU']," ",magm," ",potcar_symbols)

        """
        Write files
        """
#       write_vasp(file_path + 'slab_100/POSCAR', slab_sorted)
#       INCAR.write_file(file_path + 'slab_100/INCAR')
#       KPOINTS.write_file(file_path + 'slab_100/KPOINTS')
#       POTCAR.write_file(file_path + 'slab_100/POTCAR')
        view(slab_sorted)

        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
            elif '#PBS -q' in line:
                q_line = n
            elif '#PBS -l' in line:
                w_line = n
        PBS_N = '#PBS -N %02d_%s\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)

        if idx < 15:
            queue = 'normal'   # Change queue name if required
        else:
            queue = 'flat'
        PBS_q = '#PBS -q %s\n' % (queue)
        replace_line('jobscript_vasp.sh', q_line, PBS_q)

        walltime = '12:00:00'
        PBS_w = '#PBS -l walltime=%s\n' % (walltime)
        replace_line('jobscript_vasp.sh', w_line, PBS_w)

        destination = file_path + 'slab_100/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
#       shutil.copy(job_file, destination)    

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
     
