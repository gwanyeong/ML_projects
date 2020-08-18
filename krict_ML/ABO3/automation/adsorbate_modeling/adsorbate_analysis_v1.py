# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 10:25:52 2020

@author: gyjung
"""

import os
import csv

import pandas as pd             
import numpy as np

from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd
from ase.visualize import view


from pymatgen import MPRester
from pymatgen.io.vasp import Vasprun, Oszicar 

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################

mpr = MPRester('YOUR_MPI_KEY')
mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243

###############################################################################

entries_from_list = mpr.query(criteria = {"material_id":{"$in":mpid_list}},
                    properties = ["pretty_formula","e_above_hull"])
len(entries_from_list)  # 243

entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])

###############################################################################

adsorbates = ['OOH', 'O', 'OH']

H2O_ref = -14.21953641
H2_ref = -6.77110791

OOH_corr = 2*H2O_ref - 1.5*H2_ref
O_corr = H2O_ref - H2_ref
OH_corr = H2O_ref - 0.5*H2_ref

ads_corr = dict(zip(adsorbates,[OOH_corr, O_corr, OH_corr]))

###############################################################################

formula_list = []
for k in range(len(entries)):
    formula_list.append(entries[k]['pretty_formula'])

createFolder('results')

properties = ['formula',
              'convg_bare','tot_E_bare','DE_bare','bare_mag',
              'convg_OOH','tot_E_OOH','DE_OOH','OOH_mag','OOH_dmag',
              'convg_O', 'tot_E_O','DE_O', 'O_mag','O_dmag',
              'convg_OH','tot_E_OH','DE_OH','OH_mag','OH_dmag']

###############################################################################

with open('results/model_summary_adsorbate', 'w') as f:

    # insert the number of index for analysis here
    num = 10

    df = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                      index = [m for m in range(1,num+1)],
                      columns = properties)

    for idx in range(num):
        formula = formula_list[idx]
        df.formula[idx + 1] = formula
        createFolder('results/%03d_%s' % (idx + 1.0, formula))
        
        try:
            if os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface_np/cont/' % (idx + 1.0, formula)):
                file_path = ('%03d_%s/2nd/surface_np/cont/') % (idx + 1.0, formula)
            elif os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface/2nd/' % (idx + 1.0, formula)):
                file_path = ('%03d_%s/2nd/surface/2nd/') % (idx + 1.0, formula)
            
            bare_opt = read_vasp(file_path + 'CONTCAR')
            v_bare = Vasprun(file_path + 'vasprun.xml')
            E_bare = v_bare.ionic_steps[-1]['e_wo_entrp']

            f.writelines('%03d_%s_bare\n' % (idx + 1.0, formula))
            f.writelines('\tElectronic & ionic converged?: %s\n' % v_bare.converged) 
            f.writelines('\tFinal_energy(eV): %4.6f\n' % E_bare)

            df.convg_bare[idx + 1] = v_bare.converged
            df.tot_E_bare[idx + 1] = float('%4.4f' % E_bare)
            df.DE_bare[idx + 1] = 0
            
            oszicar = Oszicar(file_path + 'OSZICAR')
            mag_bare = oszicar.ionic_steps[-1]['mag']
            f.writelines('\tmagnetization: %4.3f\n' % mag_bare)
            df.bare_mag[idx + 1] = float('%4.3f' % mag_bare)
        
            write_xsd('results/%03d_%s/bare.xsd' % (idx + 1.0, formula), bare_opt)
          # view(bare_opt)
          
            for k in range(len(adsorbates)):  
                try:
                    slab_ini = read_vasp(file_path + '%s/POSCAR' % (adsorbates[k]))
                    slab_opt = read_vasp(file_path + '%s/CONTCAR' % (adsorbates[k]))
                    v_slab = Vasprun(file_path + '%s/vasprun.xml' % (adsorbates[k]))
        
                    E_slab = v_slab.ionic_steps[-1]['e_wo_entrp']
                    DE_ads = E_slab - E_bare - ads_corr[adsorbates[k]]                    
                                          
                    f.writelines('%03d_%s_%s\n' % (idx + 1.0, formula, adsorbates[k]))
                    f.writelines('\tElectronic & ionic converged?: %s\n' % v_slab.converged) 
                    f.writelines('\tFinal_energy(eV): %4.6f\n' % E_slab)
                    f.writelines('\tDE_ads(eV): %4.6f\n' % DE_ads)
                    
                    df.iloc[idx][5*k + 5] = v_slab.converged
                    df.iloc[idx][5*k + 6] = float('%4.4f' % E_slab)
                    df.iloc[idx][5*k + 7] = float('%4.3f' % DE_ads)
                    
                    oszicar = Oszicar(file_path + '%s/OSZICAR' % (adsorbates[k]))
                    
                    mag_ads = oszicar.ionic_steps[-1]['mag']
                    dmag = mag_ads - mag_bare
                    
                    f.writelines('\tmagnetization: %4.3f\n' % mag_ads)
                    df.iloc[idx][5*k + 8] = float('%4.3f' % mag_ads)
                    df.iloc[idx][5*k + 9] = float('%4.3f' % dmag)
                    
#                   view(slab_ini)
#                   view(slab_opt)
                    write_xsd('results/%03d_%s/%s_opt.xsd' % (idx + 1.0, formula, adsorbates[k]), slab_opt)        
                  
                except:
                    f.writelines('%03d_%s_%s files are not found !!\n' % (idx + 1.0, formula, adsorbates[k]))           
                  
        except:
            f.writelines('%03d_%s files are not found !!\n' % (idx + 1.0, formula))
            continue    

    df.to_csv('results/df_summary.csv')
