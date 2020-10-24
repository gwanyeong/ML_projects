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

from dotenv import load_dotenv

from pymatgen import MPRester
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
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

#############################################################################

# Pseudopotential choice
PP_dict = {'Sc':'Sc_sv', 'Y':'Y_sv', 'Ti':'Ti_pv', 'Zr':'Zr_sv', 'Hf':'Hf_pv',
           'V':'V_sv', 'Nb':'Nb_pv', 'Ta':'Ta_pv', 'Cr':'Cr_pv', 'Mo':'Mo_pv',
           'W':'W_sv', 'Mn':'Mn_pv', 'Tc':'Tc_pv', 'Re':'Re_pv', 'Fe':'Fe_pv',
           'Co':'Co', 'Ni':'Ni_pv', 'Cu':'Cu_pv', 'Zn':'Zn', 'Ru':'Ru_pv',
           'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'Hg':'Hg', 'Au':'Au',
           'Ir':'Ir', 'Pt':'Pt', 'Os':'Os_pv',
           'H':'H', 'O':'O'}  

#############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
                'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())

#############################################################################


INCAR = Incar.from_file('INCAR_bulk')
KPOINTS = Kpoints.from_file('KPOINTS_bulk')

with open('bulk_modeling_new.log', 'w') as f:
    start_time = time.time()
    
    num_ini = df_entries.index[0]
    num_fin = df_entries.index[-1]

    for idx in range(num_ini, num_fin):
        mp_id = df_entries['material_id'][idx]
        formula = df_entries['pretty_formula'][idx]
        struc = df_entries['structure'][idx]
        INCAR['SYSTEM'] = formula
        
        file_path  = '%03d_%s/' % (idx + 1.0, formula)
        createFolder(file_path)
     
        # getting conventional unit cell
        sga = SpacegroupAnalyzer(struc, symprec = 0.1)
        conv_struc = sga.get_conventional_standard_structure()
        conv_struc.to(filename = file_path + 'POSCAR')

        # MAGMOM setting using Exception cases
        INCAR['MAGMOM'] = [5.0, 5.0, 0.6, 0.6, 0.6]
        
        if str(conv_struc.species[0]) in alkali_elements:
            INCAR['MAGMOM'][0] = 0.6
        elif str(conv_struc.species[1]) in alkali_elements:
            INCAR['MAGMOM'][1] = 0.6

        # Hubbard U setting
        INCAR['LDAUL'] = [-1, -1, -1]
        INCAR['LDAUU'] = [0.0, 0.0, 0.0]
        INCAR['LDAUJ'] = [0.0, 0.0, 0.0]
        
        elements = [str(element) for element in conv_struc.species][:-2]
        
        for index, element in enumerate(elements):
            if element in U_elements:
                INCAR['LDAUL'][index] = 2
                INCAR['LDAUU'][index] = U_dict[element]
        
        INCAR.write_file('%03d_%s/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%03d_%s/KPOINTS' % (idx + 1.0, formula))       
        
        # Potcar setup 
        mp_calc = mpr.get_entry_by_material_id({'material_id': mp_id})
        mp_potcar_symbols = mp_calc.parameters['potcar_symbols']
        
        for i in range(len(mp_potcar_symbols)):
            mp_potcar_symbols[i] = mp_potcar_symbols[i].replace("PBE ","")
        if 'W_pv' in mp_potcar_symbols:
            mp_potcar_symbols[1] = 'W_sv'

        POTCAR = Potcar(mp_potcar_symbols)
        POTCAR.write_file(file_path + 'POTCAR')
        
        print('%02d' % (idx + 1.0)," ",formula," ",conv_struc.formula," ",len(conv_struc), " %4.3f %4.3f %4.3f" % (conv_struc.lattice.abc),
              " ",mp_potcar_symbols," ",INCAR['LDAUU']," %s" % INCAR['MAGMOM'])
        f.writelines(['%02d' % (idx + 1.0)," ",formula," ",conv_struc.formula," ",str(len(conv_struc)),
                      " %4.3f %4.3f %4.3f" % (conv_struc.lattice.abc)," %s" % mp_potcar_symbols,
                      " %s" % INCAR['LDAUU']," %s\n" % INCAR['MAGMOM']])
  
        # jobscript copy
        for n,line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
            elif '#PBS -q' in line:
                q_line = n
            
        PBS_N = '#PBS -N %03d_%s\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)
        
        if (idx - num_ini) < 50:
            queue = 'normal' # Change queue name if required
        else:
            queue = 'flat'
        PBS_q = '#PBS -q %s\n' % (queue)
        replace_line('jobscript_vasp.sh', q_line, PBS_q)
        
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, file_path)

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))


