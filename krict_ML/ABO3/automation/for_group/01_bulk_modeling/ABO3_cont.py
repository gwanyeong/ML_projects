# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 22:43:56 2020

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

from pymatgen import MPRester
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

###########################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###########################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    
############################################################################
load_dotenv(".env")
MATERIAL_API_KEY = os.getenv("MATERIAL_API_KEY")
mpr = MPRester(MATERIAL_API_KEY)

df_entries_ori = pd.read_csv('mpid_list_v2.csv')

mpid_list = df_entries_ori['mp_id'].tolist()

##############################################################################

alkali_elements = ['Sr','Ba','K','Na','Li','Rb','Cs','Be','Mg','Ca','Ba','Si']

TM_elements = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
               "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
               "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]

La_elements = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
               "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
               "Fm","Md","No","Lr"]

##############################################################################

entries = mpr.query(criteria = {"elements":{"$all":["O"], "$in":TM_elements,
                                "$nin":La_elements}, 
                                "anonymous_formula":{"A":1, "B":1, "C":3}, "nelements":3, "spacegroup.number":221,
                                "crystal_system":"cubic","spacegroup.symbol":"Pm-3m"},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])

len(entries) # 361
##############################################################################

entries = sorted(entries, key = lambda e: e['e_above_hull'])

# Remove duplicates
df_entries = pd.DataFrame(entries)
df_entries = df_entries.drop_duplicates(['pretty_formula'], keep = 'first')

df_entries.index = range(len(df_entries))
df_entries # 351 entries

# Remove the entries not belonging to ABO3
for idx in range(len(df_entries)):
    name = df_entries['pretty_formula'][idx]
    if not name.endswith('O3'):
        # print(name)
        index = df_entries[df_entries['pretty_formula'] == name].index
        df_entries = df_entries.drop(index)
df_entries.index = range(len(df_entries))
df_entries # 348 entries

# Remove the entries in the previous dataset
for idx in range(len(df_entries)):
    name = df_entries['pretty_formula'][idx]
    if name in df_entries_ori.formula.to_list():
     #   print(idx," ",name)
     #   print(df_entries[df_entries['pretty_formula'] == name])
        index = df_entries[df_entries['pretty_formula'] == name].index
        df_entries = df_entries.drop(index)
df_entries        

df_entries.index = range(len(df_entries_ori),len(df_entries_ori) + len(df_entries))
df_entries # 105 entries

df_entries.to_csv('df_entris.csv')

#############################################################################


with open('bulk_modeling_cont.log', 'w') as f:
    start_time = time.time()
    
    num_ini = 243   # df_entries.index[0]
    num_fin = 347   # df_entries.index[-1]

    for idx in range(num_ini, num_fin):
        mp_id = df_entries['material_id'][idx]
        formula = df_entries['pretty_formula'][idx]
 #      struc = df_entries['structure'][idx]
        
        file_path  = '%03d_%s/' % (idx + 1.0, formula)
        createFolder(file_path + '2nd')
      
        try:
            # vasprun.xml
            v = Vasprun(file_path + 'vasprun.xml')
            print(file_path,' Electronic & ionic converged?: %s' % v.converged)
            f.writelines([file_path,' Electronic & ionic converged?: %s\n' % v.converged])

            if v.converged:
                INCAR = Incar.from_file(file_path + 'INCAR')
                KPOINTS = Kpoints.from_file(file_path + 'KPOINTS')
                POSCAR = Poscar.from_file(file_path + 'CONTCAR')

                INCAR['IBRION'] = 2
                INCAR['ISPIN'] = 2
                INCAR['NSW'] = 200
                INCAR['ICHARG'] = 1
        
#               INCAR.write_file(file_path + '2nd/INCAR')
#               KPOINTS.write_file(file_path + '2nd/KPOINTS')
#               POSCAR.write_file(file_path + '2nd/POSCAR')
        
                # Potcar setup 
                POTCAR = Potcar.from_file(file_path + 'POTCAR')
#               POTCAR.write_file(file_path + '2nd/POTCAR')
        
                # jobscript copy
                for n,line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
                    if '#PBS -N' in line:
                        n_line = n
                    elif '#PBS -q' in line:
                        q_line = n
                    elif '#PBS -l' in line:
                        w_line = n
            
                PBS_N = '#PBS -N %03d_%s\n' % (idx + 1.0, formula)
                replace_line('jobscript_vasp.sh', n_line, PBS_N)
        
                if (idx - num_ini) < 50:
                    queue = 'normal' # Change queue name if required
                else:
                    queue = 'flat'
                PBS_q = '#PBS -q %s\n' % (queue)
                replace_line('jobscript_vasp.sh', q_line, PBS_q)

                walltime = '01:00:00'
                PBS_w = '#PBS -l walltime=%s\n' % (walltime)
                replace_line('jobscript_vasp.sh', w_line, PBS_w)
        
                destination = file_path + '2nd/'
                job_file = os.getcwd() + '/jobscript_vasp.sh'
                shutil.copy(job_file, destination)

                # CHGCAR copy
                CHGCAR = file_path + 'CHGCAR'
                if os.path.exists(file_path + '2nd/CHGCAR'):
                    print(file_path, 'CHGCAR file already exist!')
                    f.writelines([file_path, 'CHGCAR file already exist!\n'])
                else:
                    shutil.copy(CHGCAR, destination)

        except:
            print(file_path + 'error!')
            f.writelines(file_path + 'error!\n')

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))


