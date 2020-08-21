# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 22:43:56 2020

@author: gyjung
"""

import os
import shutil
import csv
import time

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

############################################################################
load_dotenv(".env")
MATERIAL_API_KEY = os.getenv("MATERIAL_API_KEY")
mpr = MPRester(MATERIAL_API_KEY)

mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243
#############################################################################
entries_from_list = mpr.query(criteria = {"material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])
len(entries_from_list)  # 243

alkali_elements = ['Sr', 'Ba', 'K', 'Na', 'Li', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Ba', 'Si']

sorted_entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])

##############################################################################

"""
entries = mpr.query(criteria = {"elements":{"$all":["O"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                            "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
                                "$nin":["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                                                   "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
                                                   "Fm","Md","No","Lr"]}, 
                    "anonymous_formula":{"A":1, "B":1, "C":3}, "nelements":3, "spacegroup.number":221,
                    "crystal_system":"cubic","spacegroup.symbol":"Pm-3m"},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])

len(entries) # 361

for i in range(len(entries)):
    print(entries[i]['pretty_formula'], end = '\t')
"""

#############################################################################

INCAR = Incar.from_file('INCAR_bulk')
KPOINTS = Kpoints.from_file('KPOINTS_bulk')

with open('bulk_modeling.log', 'w') as f:
    start_time = time.time()
    for i in range(len(sorted_entries)):
        mp_id = sorted_entries[i]['material_id']
        formula = sorted_entries[i]['pretty_formula']
        struc = sorted_entries[i]['structure']
    
        INCAR['SYSTEM'] = formula
        createFolder('%03d_%s' % (i + 1.0, formula))
     
        # getting conventional unit cell
        sga = SpacegroupAnalyzer(struc)
        conv_struc = sga.get_conventional_standard_structure()
        conv_struc.to(filename = "%03d_%s/POSCAR" % (i + 1.0, formula))
    
        # MP_incar load
        mp_incar = sorted_entries[i]['input.incar']
    
        # MAGMOM setting using Exception cases
        alkali_element = [x for x in alkali_elements if x in formula]
        if alkali_element is not None:
            INCAR['MAGMOM'] = [0.6, 5.0, 0.6, 0.6, 0.6]
        else:
            INCAR['MAGMOM'] = [5.0, 5.0, 0.6, 0.6, 0.6]
                   
        # Hubbard U setting
        try:
            INCAR['LDAU'] = mp_incar['LDAU']
            INCAR['LDAUTYPE'] = mp_incar['LDAUTYPE']
            INCAR['LDAUL'] = mp_incar['LDAUL']
            INCAR['LDAUU'] = mp_incar['LDAUU']
            INCAR['LDAUJ'] = mp_incar['LDAUJ']
        except:
            INCAR['LDAU'] = False
            INCAR['LDAUTYPE'] = 2
            INCAR['LDAUL'] = [-1, -1, -1]
            INCAR['LDAUU'] = [0.0, 0.0, 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0, 0.0]

        # Exceptional case for W        
        if 'W' in formula:
         #  print(i+1," ",formula)
            INCAR['LDAUL'] = [-1, 2, -1]
            INCAR['LDAUU'] = [0.0, 7.17, 0.0]
            INCAR['LDAUJ'] = [0.0, 0.0, 0.0]
    
        INCAR.write_file('%03d_%s/INCAR' % (i + 1.0, formula))
        KPOINTS.write_file('%03d_%s/KPOINTS' % (i + 1.0, formula))       
    
    # Potcar setup
        mp_calc = mpr.get_entry_by_material_id({'material_id': mp_id})
        mp_potcar_symbols = mp_calc.parameters['potcar_symbols']
        for idx in range(len(mp_potcar_symbols)):
            mp_potcar_symbols[idx] = mp_potcar_symbols[idx].replace("PBE ","")
        if 'W_pv' in mp_potcar_symbols:
            mp_potcar_symbols[1] = 'W_sv'
    
        # ordering for exception_list
        if mp_potcar_symbols[-1] != 'O':
     #      print(i," ", formula," ",mp_potcar_symbols)
            mp_potcar_symbols.reverse()   

        POTCAR = Potcar(mp_potcar_symbols)
        print('%d\t%s\t%s\t%s' % (i+1, formula, mp_potcar_symbols, INCAR['LDAUU']))
        f.writelines('%d\t%s\t%s\t%s\n' % (i+1, formula, mp_potcar_symbols, INCAR['LDAUU']))
        POTCAR.write_file('%03d_%s/POTCAR' % (i + 1.0, formula))
  
        # jobscript copy
        destination = '%03d_%s/' % (i + 1.0, formula)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
   

