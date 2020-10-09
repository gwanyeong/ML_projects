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
import pandas as pd
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
TM_elements = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
               "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
               "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]

nonmag_elements = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "B", "Si", "O"]

###############################################################################
createFolder('models')
df_entries = pd.read_csv('mpid_list_v2.csv')

# insert the number of index for modeling here
num_ini = 60
num_fin = 80

###############################################################################
with open('models/slab_modeling_%03d_%03d.log' % (num_ini + 1, num_fin), 'w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_slab')  

    for idx in range(num_ini,num_fin):
        formula = df_entries['formula'][idx]
        file_path = '%03d_%s/2nd/' % (idx + 1.0, formula)
        createFolder(file_path + 'surface')
   
        """
        Slab modeling
        """    
        bulk_opt = read_vasp(file_path + 'CONTCAR')
      
        slab = surface(bulk_opt, (1,0,0), 3, vacuum = 7.5)  
        slab_super = make_supercell(slab, [[2,0,0],[0,2,0],[0,0,1]])
    
        positions = slab_super.get_positions()
    
        layer_position = []
        for i in range(len(positions)):
#           print('%4.3f' % positions[i][2],end = '\t')
            if positions[i][2] not in layer_position:
                layer_position.append(positions[i][2])
    
        layer_position.sort()

        f.writelines('%03d_%s\n\tN_atoms(before) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

        natoms_bot = 0
        for atom in slab_super:
            if atom.position[2] < (layer_position[0] + 0.1):
                natoms_bot += 1
        
        if natoms_bot == 8: 
            del slab_super[[atom.position[2] < (layer_position[0] + 0.1) for atom in slab_super]]
            layer_position.remove(layer_position[0])
        elif natoms_bot == 12:
            del slab_super[[atom.position[2] > (layer_position[-1] - 0.1) for atom in slab_super]]
            layer_position.remove(layer_position[-1])
        else:
            f.writelines('N_atoms in bottomost layer - %d, which is wrong!!\n' % natoms_bot)

        f.writelines('    N_atoms(after) : %d' % (slab_super.get_global_number_of_atoms()))
    
        fix_id_list = []
        for atom in slab_super:
            if atom.position[2] < (layer_position[-3] + 0.1):
                fix_id_list.append(atom.index)
    
        f.writelines('    N_atoms(fixed) : %d\n' % (len(fix_id_list)))

        slab_super.set_constraint(FixAtoms(indices = fix_id_list))
    
        for atom in slab_super:
            atom.position[2] -= layer_position[0]

        slab_sorted = sort(slab_super)      
        view(slab_sorted)

        """
        INCAR
        """
        INCAR = Incar.from_file(file_path + 'INCAR')
        
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = 0
        INCAR['IDIPOL'] = 3
        
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
            if e in nonmag_elements:
                mag_list.append(0.6)
            elif e in TM_elements:
                mag_list.append(5.0)
            else:
                mag_list.append(5.0)

        mag_dic = dict(zip(elements, mag_list))

        magmom = []
        for el in slab_sorted.get_chemical_symbols():
            magmom.append(mag_dic[el])

#       for k in range(len(slab_sorted)):
#           print(slab_sorted[k].symbol," ",magmom[i])
            
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

        """
        Write files
        """
        write_vasp(file_path + 'surface/POSCAR', slab_sorted)
        INCAR.write_file(file_path + 'surface/INCAR')
        KPOINTS.write_file(file_path + 'surface/KPOINTS')
        POTCAR.write_file(file_path + 'surface/POTCAR')

        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
            elif '#PBS -q' in line:
                q_line = n
        PBS_N = '#PBS -N %03d_%s\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)

        queue = 'flat'   # Change queue name if required
        PBS_q = '#PBS -q %s\n' % (queue)
        replace_line('jobscript_vasp.sh', q_line, PBS_q)

        destination = file_path + 'surface/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)    

        f.writelines(['#'*80,'\n'])

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
