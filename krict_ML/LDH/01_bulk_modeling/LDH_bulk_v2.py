# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:04:16 2020

@author: gyjung
"""

import os
import shutil
import time
import fileinput
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from dotenv import load_dotenv
from pymatgen import MPRester, Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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
load_dotenv(".env")
MATERIALS_API_KEY = os.getenv("MATERIALS_API_KEY")
mpr = MPRester(MATERIALS_API_KEY)
                    
#############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
        
##############################################################################

TM_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
           'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

LDH_list = []
for TM in TM_list:
    LDH = TM + "H2O2"
    LDH_list.append(LDH)

entries = mpr.query(criteria = {"elements":{"$all":["O","H"],"$in":TM_list},
                    "anonymous_formula":{"A":1, "B":2,"C":2}, "nelements":3, "spacegroup.number":164,
                    "crystal_system":"trigonal","spacegroup.symbol":"P-3m1"},
                    properties = ["pretty_formula","structure","e_above_hull"])

entries = sorted(entries, key = lambda e: ['e_above_hull'])

df_entries = pd.DataFrame(entries)
df_entries['pretty_formula'] = [name.replace("(HO)2","H2O2") for name in df_entries['pretty_formula']]
df_entries = df_entries.drop_duplicates(['pretty_formula'], keep = 'first')
print('%d species are found among %d entries' % (len(df_entries),len(LDH_list)))

MP_list = df_entries['pretty_formula'].tolist() # ["MnH2O2","FeH2O2","CoH2O2","NiH2O2","ZnH2O2","CdH2O2"]
MP_list

MP_dict = {key:value for value, key in enumerate(MP_list)}

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

for idx, formula in enumerate(LDH_list):
    createFolder('%02d_%s' % (idx + 1.0,formula))
    
    TM=TM_list[idx]
    
    # POSCAR import/setup
    if formula in MP_list:
        index = MP_dict[formula]
        struc = df_entries['structure'][index]

        sga = SpacegroupAnalyzer(struc,0.1)
        conv_struc = sga.get_conventional_standard_structure()
        poscar = Poscar(conv_struc)
                   
    else:
        index = MP_dict['MnH2O2']
        struc_ori = df_entries['structure'][index]
        sga = SpacegroupAnalyzer(struc_ori,0.1)
        conv_struc = sga.get_conventional_standard_structure()
        conv_struc.replace(0,TM)
        poscar = Poscar(conv_struc)
        
    poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
    poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
    write_xsd('models/poscar_%02d_%s.xsd' % (idx +  1.0, formula), poscar_ase)                
    
    # INCAR setting
    INCAR['SYSTEM'] = formula
    INCAR['MAGMOM'] = [5.0, 0.6, 0.6, 0.6, 0.6]
    INCAR['ISPIN'] = 1
    INCAR['IBRION'] = -1
    INCAR['NSW'] = 0
        
    if TM in U_elements:
        INCAR['LDAU'] = '.TRUE.'
        INCAR['LDAUL'] = [2, -1, -1]
        INCAR['LDAUU'] = [U_dict[TM], 0.0, 0.0]
        INCAR['LDAUJ'] = [0.0, 0.0, 0.0]    

    else:
        INCAR['LDAU'] = '.FALSE.'
        INCAR['LDAUL'] = [-1, -1, -1]
        INCAR['LDAUU'] = [0.0, 0.0, 0.0]
        INCAR['LDAUJ'] = [0.0, 0.0, 0.0]
        
    INCAR.write_file('%02d_%s/INCAR' % (idx + 1.0, formula))
    KPOINTS.write_file('%02d_%s/KPOINTS' % (idx + 1.0, formula))
    
    # POTCAR setup
  #  POTCAR = Potcar([PP_dict[TM], PP_dict['H'], PP_dict['O']])
  #  POTCAR.write_file('%02d_%s/POTCAR' % (idx + 1.0, formula))

    print('%02d' % (idx + 1.0)," ",formula," ",conv_struc.species," ", len(conv_struc), " %4.3f %4.3f %4.3f" % (conv_struc.lattice.abc))
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
        
