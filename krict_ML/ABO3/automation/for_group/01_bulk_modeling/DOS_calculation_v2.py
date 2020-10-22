# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 11:40:47 2020

@author: gyjung
"""

import os
import shutil
import csv
import time
import fileinput
import pandas as pd

from dotenv import load_dotenv
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

df_entries = pd.read_csv('mpid_list_v2.csv')

###############################################################################

createFolder('models')

with open('models/DOS_calculation.log', 'w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_DOS')  

    num = len(df_entries)

    for idx in range(num):
        formula = df_entries['formula'][idx]
        file_path = '%03d_%s/2nd/' % (idx + 1.0, formula)
        createFolder(file_path + 'DOS')
    
        bulk_opt = read_vasp(file_path + 'CONTCAR')
        INCAR = Incar.from_file(file_path + 'INCAR')
        
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = -1
        INCAR['IBRION'] = -1
        INCAR['ICHARG'] = 11
        INCAR['NSW'] = 0
       
        POTCAR = Potcar.from_file(file_path + 'POTCAR')

        write_vasp(file_path + 'DOS/POSCAR', bulk_opt)
        INCAR.write_file(file_path + 'DOS/INCAR')
        KPOINTS.write_file(file_path + 'DOS/KPOINTS')     
        POTCAR.write_file(file_path + 'DOS/POTCAR')

        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:       
                n_line = n

        PBS_N = '#PBS -N %03d_%s_DOS\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)

        destination = file_path + 'DOS/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)    

        # CHGCAR copy
        CHGCAR = file_path + 'CHGCAR'
        if os.path.exists(file_path + 'DOS/CHGCAR'):
            print('%03d_%s DOS files already exists!' % (idx + 1.0, formula))
            f.writelines(['%03d_%s DOS files already exists!\n' % (idx + 1.0, formula)])
        else:
            shutil.copy(CHGCAR, destination)
            print('%03d_%s DOS files are generated!' % (idx + 1.0, formula))
            f.writelines(['%03d_%s DOS files are generated!' % (idx + 1.0, formula),'\n'])

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
