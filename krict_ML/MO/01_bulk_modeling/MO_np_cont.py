# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 18:19:00 2020

@author: gyjung
"""


import os
import shutil
import time
import fileinput

from dotenv import load_dotenv
from pymatgen import MPRester

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun

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

entries = mpr.query(criteria = {'elements':{'$all':["O"],
                                            "$in":['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                                                   'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                                                   'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']},
                                'anonymous_formula':{'A':1,'B':1},'nelements':2,'spacegroup.number':225,
                                'crystal_system':'cubic',
                                },
                    properties=['material_id', 'pretty_formula', 'formation_energy_per_atom',
                                'energy', 'energy_per_atom', 'structure', 'band_gap', 'input.incar',
                                'magnetic_type','total_magnetization','e_above_hull','volume','theoretical','initial_structure'])
    
print('%d entries are found!' % len(entries))

###############################################################################

MO_list = ['ScO', 'TiO', 'VO', 'CrO', 'MnO', 'FeO', 'CoO', 'NiO', 'CuO', 'ZnO',
           'YO', 'ZrO', 'NbO', 'MoO', 'TcO', 'RuO', 'RhO', 'PdO', 'AgO', 'CdO',
           'HfO', 'TaO', 'WO', 'ReO', 'OsO', 'IrO', 'PtO', 'AuO', 'HgO']

TM_list = []
for name in MO_list:
    TM_name = name[:-1]
    TM_list.append(TM_name)
    
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
for idx, formula in enumerate(MO_list):
    for k in range(len(entries)):
        if formula == entries[k]['pretty_formula']:
            abc = entries[k]['structure'].lattice.abc
            print(idx," ", k," ",formula," ","%4.3f %4.3f %4.3f " % (abc))
"""

with open('np_cont_modeling.log', 'w') as f:
    start_time = time.time()

    for idx, formula in enumerate(MO_list):
        
        file_path = '%02d_%s/' % (idx + 1.0, formula)
        createFolder(file_path + 'cont')
        
        
        try:
            INCAR = Incar.from_file(file_path + 'INCAR')
            KPOINTS = Kpoints.from_file(file_path + 'KPOINTS')
            
            POSCAR = Poscar.from_file(file_path + 'CONTCAR')
            POSCAR.write_file(file_path + 'cont/POSCAR')
            
            INCAR['IBRION'] = 2
            INCAR['ISPIN'] = 2
            INCAR['MAGMOM'] = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6]
            INCAR['NSW'] = 200
            INCAR['ICHARG'] = 1
            
            INCAR.write_file(file_path + 'cont/INCAR')
            KPOINTS.write_file(file_path + 'cont/KPOINTS')
            
            # Potcar setup
            POTCAR = Potcar.from_file(file_path + 'POTCAR')
            POTCAR.write_file(file_path + 'cont/POTCAR')
            
            # vasprun.xml
            v = Vasprun(file_path + 'vasprun.xml')
            print(file_path,' Electronic & ionic converged?: %s' % v.converged)
            f.writelines([file_path, ' Electronic & ionic converged?: %s\n' % v.converged])
            
            # jobscript copy
            for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
                if '#PBS -N' in line:
                    n_line = n
            PBS_N = '#PBS -N %02d_%s\n' % (idx + 1.0, formula)
            replace_line('jobscript_vasp.sh', n_line, PBS_N)
            
            destination = file_path + 'cont/'
            job_file = os.getcwd() + '/jobscript_vasp.sh'
            shutil.copy(job_file, destination)
            
            # CHGCAR copy
            CHGCAR = file_path + 'CHGCAR'
            if os.path.exists(file_path + 'cont/CHGCAR'):
                print(file_path, 'CHGCAR file already exist!!')
                f.writelines([file_path, 'CHGCAR file already exist!!\n'])
                
            else:
                shutil.copy(CHGCAR, destination)
                
        except:
            print(file_path + 'calculations are not properly finished!')
            f.writelines(file_path + 'calculations are not properly finished!\n')
            
    end_time = time.time()
    print('Execution time for script (sec) : %5.1f\n' % (end_time - start_time))
    f.writelines('Execution time for script (sec) : %5.1f\n' % (end_time - start_time))
            