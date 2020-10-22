
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:13:25 2020
@author: gyjung
"""

import os
import shutil
import time
import fileinput
import warnings
warnings.filterwarnings('ignore')

import pandas as pd

from dotenv import load_dotenv
from pymatgen import MPRester

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ase.visualize import view
from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd

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

load_dotenv(".env")
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')
mpr = MPRester(MATERIAL_API_KEY)

###############################################################################

M2O3_list = ['Sc2O3', 'Ti2O3', 'V2O3', 'Cr2O3', 'Mn2O3', 'Fe2O3', 'Co2O3', 'Ni2O3', 'Cu2O3', 'Zn2O3',
           'Y2O3', 'Zr2O3', 'Nb2O3', 'Mo2O3', 'Tc2O3', 'Ru2O3', 'Rh2O3', 'Pd2O3', 'Ag2O3', 'Cd2O3',
           'Hf2O3', 'Ta2O3', 'W2O3', 'Re2O3', 'Os2O3', 'Ir2O3', 'Pt2O3', 'Au2O3', 'Hg2O3']

TM_list = []
for name in M2O3_list:
    TM_name = name.replace("2O3","")
    TM_list.append(TM_name)

###############################################################################

properties = ['pretty_formula','structure','e_above_hull']

entries = mpr.query(criteria = {'elements':{'$all':["O"],
                                            "$in":['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                                                   'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                                                   'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']},
                                'anonymous_formula':{'A':2,'B':3},'nelements':2,'spacegroup.number':167,
                                'crystal_system':'trigonal','spacegroup.symbol':'R-3c'
                                },
                    properties=properties)
    
print('%d species are found among %d entries' % (len(entries),len(TM_list)))

entries = sorted(entries, key = lambda e: e['e_above_hull'])


for entry in entries:
    # print(entry['pretty_formula'])
    if entry['pretty_formula'] not in M2O3_list:
        entries.remove(entry)

df = pd.DataFrame(entries)
df = df.drop_duplicates(['pretty_formula'], keep = 'first')

###############################################################################

# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2, 'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())

# Pseudopotential choice

PP_dict = {'Sc':'Sc_sv', 'Y':'Y_sv', 'Ti':'Ti_pv', 'Zr':'Zr_sv', 'Hf':'Hf_pv',
           'V':'V_sv', 'Nb':'Nb_pv', 'Ta':'Ta_pv', 'Cr':'Cr_pv', 'Mo':'Mo_pv',
           'W':'W_sv', 'Mn':'Mn_pv', 'Tc':'Tc_pv', 'Re':'Re_pv', 'Fe':'Fe_pv',
           'Co':'Co', 'Ni':'Ni_pv', 'Cu':'Cu_pv', 'Zn':'Zn', 'Ru':'Ru_pv',
           'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'Hg':'Hg', 'Au':'Au',
           'Ir':'Ir', 'Pt':'Pt', 'Os':'Os_pv',
           'H':'H', 'O':'O'}    


createFolder('models')

"""
for idx, formula in enumerate(M2O3_list):
    for k in range(len(entries)):
        if formula == entries[k]['pretty_formula']:
            abc = entries[k]['structure'].lattice.abc
            print(idx," ", k," ",formula," ","%4.3f %4.3f %4.3f " % (abc))
"""


with open('M2O3_bulk_modeling.log', 'w') as f:
    start_time = time.time()
    INCAR = Incar.from_file('INCAR')
    KPOINTS = Kpoints.from_file('KPOINTS')

    formula_ini = df['pretty_formula'][0]
    struc_ini = df['structure'][0] # Fe2O3
        
    f.writelines("Initial model: %s\t" % formula_ini)
    f.writelines("lattice : %4.3f, %4.3f %4.3f\n" % (struc_ini.lattice.abc))

    for idx, formula in enumerate(M2O3_list):
        createFolder('%02d_%s' % (idx + 1.0, formula))

        struc_ini.replace(0, TM_list[idx])
        struc_ini.replace(1, TM_list[idx])
        struc_ini.replace(2, TM_list[idx])
        struc_ini.replace(3, TM_list[idx])
         
        # getting conventional unit cell
        sga = SpacegroupAnalyzer(struc_ini)
        conv_struc = sga.get_conventional_standard_structure()
        
        poscar = Poscar(conv_struc)
        poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
    
        poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
        write_xsd('models/poscar_%02d_%s.xsd' % (idx + 1.0, formula), poscar_ase)
        
        for k in range(len(df)):
            if formula == df['pretty_formula'][k]:
            
                struc = df['structure'][k]
                sga = SpacegroupAnalyzer(struc)
                # print(idx+1.0," ",formula," ",len(sga.get_conventional_standard_structure()))
                
                if len(sga.get_conventional_standard_structure()) != 30:
                    print('Error for conv. struc. of %s' % formula)
                    f.writelines('Error for conv. struc. of %s\n' % formula)
                    break            
                conv_struc = sga.get_conventional_standard_structure()
                poscar = Poscar(conv_struc)
                poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
                poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
   #            view(poscar_ase)
                write_xsd('models/poscar_%02d_%s.xsd' % (idx + 1.0, formula), poscar_ase)
    
        print('%02d' % (idx + 1.0)," ",formula," ",poscar.structure.formula," ",len(conv_struc), " %4.3f %4.3f %4.3f" % (conv_struc.lattice.abc))
        f.writelines(['%02d' % (idx + 1.0)," ",formula," ",poscar.structure.formula," ",str(len(conv_struc)), " %4.3f %4.3f %4.3f\n" % (conv_struc.lattice.abc)])

        # INCAR setup
        INCAR['SYSTEM'] = formula
      # INCAR['MAGMOM'] = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6]
        INCAR['ISPIN'] = 1
        INCAR['IBRION'] = -1
        INCAR['NSW'] = 0
    
        TM = formula.replace("2O3","")
    
        if TM in U_elements:
            INCAR['LDAU'] = '.TRUE.'
            INCAR['LDAUL'] = [2, -1]
            INCAR['LDAUU'] = [U_dict[TM], 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0]
        
        else:
            INCAR['LDAU'] = '.FALSE.'
            INCAR['LDAUL'] = [-1, -1]
            INCAR['LDAUU'] = [0.0, 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0]
        
        INCAR.write_file('%02d_%s/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%02d_%s/KPOINTS' % (idx + 1.0, formula))
    
        # Potcar setup
        POTCAR = Potcar([PP_dict[TM], PP_dict['O']])
        POTCAR.write_file('%02d_%s/POTCAR' % (idx + 1.0, formula))
      
        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
        PBS_N = '#PBS -N %02d_%s_np\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)
    
        destination = '%02d_%s/' % (idx + 1.0, formula)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)
    
    end_time = time.time()
    print('Execution time for script (sec) : %5.1f\n' % (end_time - start_time))
    f.writelines('Execution time for script (sec) : %5.1f\n' % (end_time - start_time))
