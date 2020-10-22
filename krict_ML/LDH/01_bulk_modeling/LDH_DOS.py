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
LDH_list = ["ScH2O2","TiH2O2","VH2O2", "CrH2O2","MnH2O2","FeH2O2","CoH2O2","NiH2O2","CuH2O2","ZnH2O2",
            "YH2O2","ZrH2O2","NbH2O2","MoH2O2","TcH2O2","RuH2O2","RhH2O2","PdH2O2","AgH2O2","CdH2O2",
            "HfH2O2","TaH2O2","WH2O2","ReH2O2","OsH2O2","IrH2O2","PtH2O2","AuH2O2","HgH2O2"]

###############################################################################
createFolder('models')

with open('models/DOS_calculation.log','w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_DOS')
    
    for idx, formula in enumerate(LDH_list):
        file_path = '%02d_%s/opt/' % (idx + 1.0, formula)
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
