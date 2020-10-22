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

from dotenv import load_dotenv

from ase.io.vasp import read_vasp, write_vasp

from pymatgen import MPRester
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

##############################################################################

MHO2_list = ["ScHO2","TiHO2","VHO2","CrHO2","MnHO2","FeHO2","CoHO2","NiHO2","CuHO2","ZnHO2",
            "YHO2","ZrHO2","NbHO2","MoHO2","TcHO2","RuHO2","RhHO2","PdHO2","AgHO2","CdHO2",
            "HfHO2","TaHO2","WHO2","ReHO2","OsHO2","IrHO2","PtHO2","AuHO2","HgHO2"]

TM_list=[]
for name in MHO2_list:
    TM_name = name.replace("HO2","")
    TM_list.append(TM_name)

###############################################################################

# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
                'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())


# Pseudopotential choice
PP_dict = {'Sc':'Sc_sv', 'Y':'Y_sv', 'Ti':'Ti_pv', 'Zr':'Zr_sv', 'Hf':'Hf_pv',
           'V':'V_sv', 'Nb':'Nb_pv', 'Ta':'Ta_pv', 'Cr':'Cr_pv', 'Mo':'Mo_pv',
           'W':'W_sv', 'Mn':'Mn_pv', 'Tc':'Tc_pv', 'Re':'Re_pv', 'Fe':'Fe_pv',
           'Co':'Co', 'Ni':'Ni_pv', 'Cu':'Cu_pv', 'Zn':'Zn', 'Ru':'Ru_pv',
           'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'Hg':'Hg', 'Au':'Au',
           'Ir':'Ir', 'Pt':'Pt', 'Os':'Os_pv',
           'H':'H', 'O':'O'}    

##############################################################################

createFolder('models')

with open('models/DOS_calculation.log', 'w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_DOS')  

for idx, formula in enumerate(MHO2_list):
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
       
        POTCAR = Potcar.from_file(file_path + 'POTCAR')

        write_vasp(file_path + 'DOS/POSCAR', bulk_opt)
        INCAR.write_file(file_path + 'DOS/INCAR')
        KPOINTS.write_file(file_path + 'DOS/KPOINTS')     
        POTCAR.write_file(file_path + 'DOS/POTCAR')

        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:       
                n_line = n

        PBS_N = '#PBS -N %02d_%s_DOS\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)

        destination = file_path + 'DOS/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)    

        # CHGCAR copy
        CHGCAR = file_path + 'CHGCAR'
        if os.path.exists(file_path + 'DOS/CHGCAR'):
                print(file_path, 'CHGCAR file already exist!!')
                f.writelines([file_path, 'CHGCAR file already exist!!\n'])
        else:
                shutil.copy(CHGCAR, destination)

end_time = time.time()
f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
