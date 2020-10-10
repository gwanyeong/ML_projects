# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 23:14:38 2020
@author: hjkim
"""

import os
import shutil
import csv
import time
import fileinput

from dotenv import load_dotenv
from pymatgen import MPRester
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

LDH_list = ["ScH2O2","TiH2O2","VH2O2","CrH2O2","MnH2O2","FeH2O2","CoH2O2","NiH2O2","CuH2O2","ZnH2O2",
            "YH2O2","ZrH2O2","NbH2O2","MoH2O2","TcH2O2","RuH2O2","RhH2O2","PdH2O2","AgH2O2","CdH2O2",
            "HfH2O2","TaH2O2","WH2O2","ReH2O2","OsH2O2","IrH2O2","PtH2O2","AuH2O2","HgH2O2"]

TM_list=[]
for name in LDH_list:
    TM_name = name.replace("H2O2","")
    TM_list.append(TM_name)

MP_list=["MnH2O2","FeH2O2","CoH2O2","NiH2O2","ZnH2O2","CdH2O2"]

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
    
    TM=formula.replace("H2O2","")
    
    if formula in MP_list:
        entries = mpr.query(criteria = {"elements":{"$all":["O","H"], "$in":[TM]},
                    "anonymous_formula":{"A":1, "B":2,"C":2}, "nelements":3, "spacegroup.number":164,
                    "crystal_system":"trigonal","spacegroup.symbol":"P-3m1"},
                    properties = ["material_id","pretty_formula","cifs","structure"])
        for i in range(len(entries)):
            struc=entries[i]['structure']
            
            sga=SpacegroupAnalyzer(struc,0.1)
            conv_struc=sga.get_conventional_standard_structure()
            
            poscar = Poscar(conv_struc)
            poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
            poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
            write_xsd('models/poscar_%02d_%s.xsd' % (idx +  1.0, formula), poscar_ase)
            id=entries[i]['material_id']
            formula_MP=entries[i]['pretty_formula']
            print(id,formula_MP)
            
        INCAR['SYSTEM'] = formula
        INCAR['MAGMOM'] = [5.0, 0.6, 0.6, 0.6, 0.6]
        INCAR['ISPIN'] = 1
        INCAR['IBRION'] = -1
        INCAR['NSW'] = 0
        
        if TM in U_elements:
            INCAR['LDAU'] = '.TRUE.'
            INCAR['LDAUL'] = [2, -1, -1]
            INCAR['LDAUU'] = [U_dic[TM], 0.0, 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0, 0.0]    

        else:
            INCAR['LDAU'] = '.FALSE.'
            INCAR['LDAUL'] = [-1, -1, -1]
            INCAR['LDAUU'] = [0.0, 0.0, 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0, 0.0]
        
        INCAR.write_file('%02d_%s/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%02d_%s/KPOINTS' % (idx + 1.0, formula))
    
        # Potcar setup
        POTCAR = Potcar([PP_dict[TM], PP_dict['H'], PP_dict['O']])
        POTCAR.write_file('%02d_%s/POTCAR' % (idx + 1.0, formula))
        
    else:
        poscar_ini = Poscar.from_file('POSCAR')                 
        struc_ori = poscar_ini.structure
    
        struc_ori.replace(2,TM_list[idx])
    
        poscar = Poscar(struc_ori)
        poscar.write_file('%02d_%s/POSCAR' % (idx + 1.0, formula))
        poscar_ase = read_vasp('%02d_%s/POSCAR' % (idx + 1.0, formula))
        write_xsd('models/poscar_%02d_%s.xsd' % (idx +  1.0, formula), poscar_ase)
    
        INCAR['SYSTEM'] = formula
        INCAR['MAGMOM'] = [0.6, 0.6, 5.0, 0.6, 0.6]
        INCAR['ISPIN'] = 1
        INCAR['IBRION'] = -1
        INCAR['NSW'] = 0
            
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
        
        
