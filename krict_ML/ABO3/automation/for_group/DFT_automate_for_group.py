# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 22:43:56 2020

@author: gyjung
"""


import os
# from dotenv import load_dotenv
import shutil

from pymatgen import MPRester
import csv

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


# load_dotenv(".env")
# MATERIALS_KEY = os.getenv("MATERIALS_KEY")

mpr = MPRester('mwEnJkSB9Cusgdf1')

mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243개

entries_from_list = mpr.query(criteria = {"elements":{"$all":["O"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                            "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
                                "$nin":["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                                                   "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
                                                   "Fm","Md","No","Lr"]}, 
                    "anonymous_formula":{"A":1, "B":1, "C":3}, "nelements":3, "spacegroup.number":221,
                    "crystal_system":"cubic","spacegroup.symbol":"Pm-3m",
                    "material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])
len(entries_from_list)  # 243개

entries_inc_alkali = mpr.query(criteria = {"elements":{"$all":["O"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                            "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
                                            "$in":["Sr", "Ba", "K", "Na", "Li", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Si"],    
                                "$nin":["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                                                   "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
                                                   "Fm","Md","No","Lr"]}, 
                    "anonymous_formula":{"A":1, "B":1, "C":3}, "nelements":3, "spacegroup.number":221,
                    "crystal_system":"cubic","spacegroup.symbol":"Pm-3m",
                    "material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])

len(entries_inc_alkali)

entries_alkali_list = []

for i in range(len(entries_inc_alkali)):
    entries_alkali_list.append(entries_inc_alkali[i]['pretty_formula'])

entries_inc_W = mpr.query(criteria = {"elements":{"$all":["O","W"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                            "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
                                      "$nin":["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                                                   "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
                                                   "Fm","Md","No","Lr"]}, 
                    "anonymous_formula":{"A":1, "B":1, "C":3}, "nelements":3, "spacegroup.number":221,
                    "crystal_system":"cubic","spacegroup.symbol":"Pm-3m",
                    "material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical"])

len(entries_inc_W)

entries_W_list = []

for i in range(len(entries_inc_W)):
    entries_W_list.append(entries_inc_W[i]['pretty_formula'])


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

len(entries) # 361개

for i in range(len(entries)):
    print(entries[i]['pretty_formula'], end = '\t')
"""

sorted_entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])


"""
for i in range(len(sorted_entries)):
    print(sorted_entries[i]['pretty_formula'], " %4.3f" % sorted_entries[i]['e_above_hull'])
"""


from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

INCAR = Incar.from_file('INCAR_bulk')
KPOINTS = Kpoints.from_file('KPOINTS_bulk')

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
    try:
        if formula in entries_alkali_list:
            INCAR['MAGMOM'] = [0.6, 5.0, 0.6, 0.6, 0.6]
        else:
            INCAR['MAGMOM'] = [5.0, 5.0, 0.6, 0.6, 0.6]
    except:
       pass
                   
    # Hubbard U setting
    # print(i," ", formula)
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
        
    if formula in entries_W_list:
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
    exception_list = [2, 22, 58]
    if i in exception_list:
        mp_potcar_symbols.reverse()   

    POTCAR = Potcar(mp_potcar_symbols)
    print(i+1," ",formula," ",mp_potcar_symbols," ",INCAR['LDAUU'])
    POTCAR.write_file('%03d_%s/POTCAR' % (i + 1.0, formula))
  
    # jobscript copy
    destination = '%03d_%s/' % (i + 1.0, formula)
    job_file = os.getcwd() + '/jobscript_vasp.sh'
    shutil.copy(job_file, destination)


