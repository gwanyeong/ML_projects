# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:13:09 2020

@author: gyjung
"""


# import os
import csv
from pymatgen import MPRester


from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.outputs import Oszicar


import pandas as pd

with open('form_E_ref.csv', 'r') as f:
    ref_data = pd.read_csv(f)
    ref_dict = dict.fromkeys(ref_data['element'])
    for i in range(len(ref_data)):
        ref_dict.update({ref_data['element'][i]:ref_data['E'][i]})

# add O2_corr ref. energy
ref_dict['O2_corr'] = -8.45572061

# define U correction values
U_corr_dict = {'V':1.682, 'Cr':2.013, 'Mn':1.68085, 'Fe':2.733,
               'Co':1.874, 'Ni':2.164, 'Mo':3.531, 'W':5.159}


mpr = MPRester('YOUR_MPI_KEY')
mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243개


entries_from_list = mpr.query(criteria = {"material_id":{"$in":mpid_list}},
                    properties = ["pretty_formula","e_above_hull"])
len(entries_from_list)  # 243개    


sorted_entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])


with open('summary', 'w') as f:
    f.writelines('No\tformula\ttotal E\tMagm\tF.E.(eV/f.u.)\tVolume\n')
    for i in range(0,243):
        formula = sorted_entries[i]['pretty_formula']
        
        try:
            oszicar = Oszicar('%03d_%s/2nd/OSZICAR'% (i + 1.0, formula))
            tot_E = oszicar.ionic_steps[-1]['F']
            tot_mag = oszicar.ionic_steps[-1]['mag']
        except FileNotFoundError:
            print('%03d_%s/2nd files are not found' % (i + 1.0, formula))
            f.writelines('%03d_%s/2nd files are not found\n' % (i + 1.0, formula))
            continue

        try:
            v = Vasprun('%03d_%s/2nd/vasprun.xml' % (i + 1.0, formula))
    
            volume = v.as_dict()['output']['crystal']['lattice']['volume']
            #print(volume)    
            el = v.as_dict()['elements']
            el.remove('O')
    
            fE = tot_E - ref_dict[el[0]] - ref_dict[el[1]] - 1.5*ref_dict['O2_corr']
    
            if el[0] or el[1] in U_corr_dict:
                for x in range(len(el)):
                    if el[x] in U_corr_dict:
                        fE = fE - U_corr_dict[el[x]]

            #print(i+1," ", formula," ",tot_E," ",tot_mag," ",fE," ",volume) 
            f.writelines('%d\t%s\t%4.4f\t%4.3f\t%4.3f\t%4.3f\n' % (i+1, formula, tot_E, tot_mag, fE, volume))

        except FileNotFoundError:
          #  print('%03d_%s/2nd files are not found' %(i + 1.0, formula))
          #  f.writelines('%03d_%s/2nd files are not found\n' %(i + 1.0, formula))
            continue

        
