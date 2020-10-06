# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 23:14:38 2020

@author: gyjung
"""

import os
import shutil
import time
import fileinput

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar

#from ase.visualize import view
from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd

start_time = time.time()

###########################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

#############################################################################
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

##############################################################################

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

INCAR = Incar.from_file('INCAR')
KPOINTS = Kpoints.from_file('KPOINTS')

createFolder('models')

for idx, formula in enumerate(MHO2_list):
    createFolder('%02d_%s' % (idx + 1.0,formula))
    
    poscar_ini = Poscar.from_file('POSCAR')                 
    struc_ori = poscar_ini.structure
    
    struc_ori.replace(2,TM_list[idx])
    struc_ori.replace(3,TM_list[idx])
    
    poscar = Poscar(struc_ori)
    poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
    
    poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
    write_xsd('models/poscar_%02d_%s.xsd' % (idx +  1.0, formula), poscar_ase)
    
    INCAR['SYSTEM'] = formula
   # INCAR['MAGMOM'] = [0.6, 0.6, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6]
    INCAR['ISPIN'] = 1
    INCAR['IBRION'] = -1
    INCAR['NSW'] = 0
        
    TM = formula.replace("HO2","")
    
    if TM in U_elements:
        INCAR['LDAU'] = '.TRUE.'
        INCAR['LDAUL'] = [-1, 2, -1]
        INCAR['LDAUU'] = [0.0, U_dict[TM], 0.0]
        INCAR['LDAUJ'] = [0.0, 0.0, 0.0]

    else:
        INCAR['LDAU'] = '.FALSE.'
        INCAR['LDAUL'] = [-1, -1, -1]
        INCAR['LDAUU'] = [0.0, 0.0, 0.0]
        INCAR['LDAUJ'] = [0.0, 0.0, 0.0]
        
    INCAR.write_file('%02d_%s/INCAR' % (idx + 1.0, formula))
    KPOINTS.write_file('%02d_%s/KPOINTS' % (idx + 1.0, formula))
    
    # Potcar setup
    POTCAR = Potcar([PP_dict['H'], PP_dict[TM], PP_dict['O']])
    POTCAR.write_file('%02d_%s/POTCAR' % (idx + 1.0, formula))
    
    # jobscript copy
    
    for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
        if '#PBS -N' in line:
            n_line = n
    PBS_N = '#PBS -N %02d_%s_np\n' % (idx+1.0, formula)    
    replace_line('jobscript_vasp.sh', n_line, PBS_N)
    
    destination = '%02d_%s/' % (idx + 1.0, formula)
    job_file = os.getcwd() + '/jobscript_vasp.sh'
    shutil.copy(job_file, destination)
    
end_time = time.time()
print('Execution time for script (sec) : %5.1f\n' % (end_time - start_time))
        
        
