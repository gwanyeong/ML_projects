# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 17:51:19 2020
@author: gyjung
"""

import os
import shutil
import time
import fileinput
import warnings
warnings.filterwarnings('ignore')

from dotenv import load_dotenv
from pymatgen import MPRester

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


from ase.visualize import view
from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)
        
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

load_dotenv(".env")
MATERIAL_API_KEY = os.getenv("MATERIAL_API_KEY")
mpr = MPRester(MATERIAL_API_KEY)


entries = mpr.query(criteria = {"elements":{"$all":["O"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                   "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                   "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]}, 
                                   "anonymous_formula":{"A":1, "B":2}, "nelements":2, "spacegroup.number":164,"crystal_system":"trigonal", # "spacegroup.symbol":"P42/mnm"
                                   },
                    properties = ["pretty_formula",
                                  "structure",
                                  ])

print('%d entries are found!' % len(entries))

for i in range(len(entries)):
    print(i," ",entries[i]['pretty_formula'])
# print(entries)
    
LTMO_list = ["ScO2","TiO2","VO2",
             "CrO2","MnO2","FeO2","CoO2","NiO2","CuO2","ZnO2",
             "YO2","ZrO2","NbO2","MoO2","TcO2","RuO2","RhO2","PdO2","AgO2","CdO2",
             "HfO2","TaO2","WO2","ReO2","OsO2","IrO2","PtO2","AuO2","HgO2"]

TM_list=[]
for name in LTMO_list:
    TM_name = name.replace("O2","")
    TM_list.append(TM_name)
    
    
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
"""
for entry in entries[0]['structure']:
    print(entry)

for idx, formula in enumerate(LTMO_list):
    for k in range(len(entries)):
        if formula == entries[k]['pretty_formula']:
            print(idx," ",k," ",formula)
"""        
createFolder('models')

with open('LTMO_bulk_modeling.log', 'w') as f:
    start_time = time.time()
    INCAR = Incar.from_file('INCAR')
    KPOINTS = Kpoints.from_file('KPOINTS')
    
    for idx, formula in enumerate(LTMO_list):

        file_path = '%02d_%s/' % (idx + 1.0, formula)
        createFolder(file_path)
        
        struc_ori = entries[0]['structure'] # PtO2                 
        
        struc_ori.replace(0,TM_list[idx])
                
        conv_struc = struc_ori
               
        poscar = Poscar(conv_struc)
        poscar.write_file(file_path + 'POSCAR')        
        
        poscar_ase = read_vasp(file_path + 'POSCAR')
        write_xsd('models/poscar_%02d_%s.xsd' % (idx + 1.0, formula), poscar_ase)
        
        # view(poscar_ase)
        
        for k in range(len(entries)):
            if formula == entries[k]['pretty_formula']:
                
                struc = entries[k]['structure']
                
               # getting conventional unit cell
                sga = SpacegroupAnalyzer(struc, symprec = 0.1)
                
                # print(idx," ",k," ",sga.get_conventional_standard_structure().formula)
                
                if len(sga.get_conventional_standard_structure()) != 3:
                    break
                
                conv_struc = sga.get_conventional_standard_structure()
               # print(idx," ",k," ",formula," ",conv_struc.species," ",len(struc), len(conv_struc))
                
                poscar = Poscar(conv_struc)
                poscar.write_file(file_path + 'POSCAR')           
                poscar_ase = read_vasp(file_path + 'POSCAR')
          
                if poscar_ase[0].position[2] > 0.0:
                #print(idx," ",formula)
                    z_coord = poscar_ase[0].position[2]
                    for atom in poscar_ase:
                        atom.position[2] -= z_coord
                 #   view(poscar_ase)
                    # print(poscar_ase.get_scaled_positions())
                
                write_vasp(file_path + 'POSCAR', poscar_ase)
                poscar = Poscar.from_file(file_path + 'POSCAR')
                conv_struc = poscar.structure
                write_xsd('models/poscar_%02d_%s.xsd' % (idx + 1.0, formula), poscar_ase)
                
        abc = "%4.3f  %4.3f  %4.3f" % conv_struc.lattice.abc
        print('%02d' % (idx + 1.0)," ",formula," ",conv_struc.species," ",len(conv_struc)," %s" % abc)               
        f.writelines(['%02d' % (idx + 1.0)," ",formula," ",str(conv_struc.species)," ",str(len(conv_struc))," %s\n" % abc])   
     #   view(poscar_ase)        
    
        INCAR['SYSTEM'] = formula           
#       INCAR['MAGMOM'] = [5.0, 0.6, 0.6]
        
        TM = formula.replace("O2","")
        if TM in U_elements:
            INCAR['LDAU'] = '.TRUE.'
            INCAR['LDAUL'] = [2, -1]
            INCAR['LDAUU'] = [U_dict[TM], 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0]
        
        INCAR.write_file('%02d_%s/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%02d_%s/KPOINTS' % (idx + 1.0, formula))
        
        # Potcar setup
        POTCAR = Potcar([PP_dict[TM],PP_dict['O']])
        POTCAR.write_file('%02d_%s/POTCAR' % (idx + 1.0, formula))
    
        # jobscript copy
        for n,line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
    
        PBS_N = '#PBS -N %02d_%s_np\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)
    
        destination = '%02d_%s/' % (idx + 1.0, formula)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)
    
    end_time = time.time()
    print('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))

