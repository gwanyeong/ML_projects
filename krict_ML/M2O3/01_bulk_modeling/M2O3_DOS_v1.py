#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 11:40:47 2020
@author: gyjung
"""

import os
import shutil
import time
import fileinput
import warnings
warnings.filterwarnings('ignore')

from ase.io.vasp import read_vasp, write_vasp

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

###############################################################################

M2O3_list = ["Sc2O3", "Ti2O3","V2O3", "Cr2O3","Mn2O3","Fe2O3","Co2O3","Ni2O3","Cu2O3","Zn2O3",
            "Y2O3","Zr2O3","Nb2O3","Mo2O3","Tc2O3","Ru2O3","Rh2O3","Pd2O3","Ag2O3","Cd2O3",
            "Hf2O3","Ta2O3","W2O3","Re2O3","Os2O3","Ir2O3","Pt2O3","Au2O3","Hg2O3"]


###############################################################################
createFolder('models')

with open('models/DOS_calculation.log','w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_DOS')
    
    for idx, formula in enumerate(M2O3_list):
        file_path = '%02d_%s/cont/' % (idx + 1.0, formula)
        createFolder(file_path + 'DOS')
        
        bulk_opt = read_vasp(file_path + 'CONTCAR')
        INCAR = Incar.from_file(file_path + 'INCAR')
        
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = -1
        INCAR['IBRION'] = -1
        INCAR['ICHARG'] = 11
        INCAR['NSW'] = 0
        INCAR['IDIPOL'] = 0
        
        POTCAR = Potcar.from_file(file_path + 'POTCAR')
        
        write_vasp(file_path + 'DOS/POSCAR',bulk_opt)
        INCAR.write_file(file_path + 'DOS/INCAR')
        KPOINTS.write_file(file_path + 'DOS/KPOINTS')
        POTCAR.write_file(file_path + 'DOS/POTCAR')
        
        #jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
            elif '#PBS -l' in line:
                w_line = n
        
        PBS_N = '#PBS -N %02d_%s_DOS\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)
        
        walltime = '06:00:00'
        PBS_w = '#PBS -l walltime=%s\n' % (walltime)
        replace_line('jobscript_vasp.sh', w_line, PBS_w)
        
        destination = file_path + 'DOS/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)
        
        #CHGCAR copy
        CHGCAR = file_path + 'CHGCAR'
        if os.path.exists(file_path + 'DOS/CHGCAR'):
            print('%02d_%s DOS files already exists!' % (idx + 1.0, formula))
            f.writelines(['%02d_%s DOS files already exists!\n' % (idx + 1.0, formula), '\n'])
        else:
            shutil.copy(CHGCAR, destination)
            print('%02d_%s DOS files are generated!' % (idx + 1.0, formula))
            f.writelines(['%02d_%s DOS files are generated!' % (idx + 1.0, formula), '\n'])
            
    end_time = time.time()
    f.writelines('Excution time for script (sec) : %6.1f\n' % (end_time - start_time))

