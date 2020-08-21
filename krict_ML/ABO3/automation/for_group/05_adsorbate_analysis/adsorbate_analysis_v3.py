# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 10:25:52 2020

@author: gyjung
"""

import os
import csv
import pickle
import time

import pandas as pd             
import numpy as np

from dotenv import load_dotenv

from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd
from ase.visualize import view


from pymatgen import MPRester
from pymatgen.io.vasp import Vasprun, Oszicar, Outcar

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
load_dotenv('.env')
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')

mpr = MPRester(MATERIAL_API_KEY)

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
              'convg_bare','tot_E_bare','DE_bare','bare_mag','nsteps_bare','time_bare',
              'convg_OOH','tot_E_OOH','DE_OOH','OOH_mag','OOH_dmag','nsteps_OOH','time_OOH',
              'convg_O', 'tot_E_O','DE_O', 'O_mag','O_dmag','nsteps_O','time_O',
              'convg_OH','tot_E_OH','DE_OH','OH_mag','OH_dmag','nsteps_OH','time_OH']

magm_properties = ['magm_bare', 'magm_OOH', 'magm_O', 'magm_OH']

###############################################################################


# insert the number of index for analysis here
num_ini = 10
num_fin = 20
num = num_fin - num_ini

with open('results/results_summary_%03d_%03d' % (num_ini+1, num_fin), 'w') as f:     # put AFM, def in the path if required
    start_time = time.time()
   
    df = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                      index = [m for m in range(num_ini+1 ,num_fin+1)],
                      columns = properties)
    df_mag = pd.DataFrame(np.array([['NaN' for p in range(len(magm_properties))] for q in range(num)]),
                          index = [r for r in range(num_ini+1, num_fin+1)],
                          columns = magm_properties)

    for idx in range(num_ini, num_fin):
        formula = formula_list[idx]
        df.formula[idx + 1] = formula
        createFolder('results/%03d_%s' % (idx + 1.0, formula))  # put AFM, def in the path if required
        
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
            df.nsteps_bare[idx + 1] = len(v_bare.ionic_steps)
            
            oszicar = Oszicar(file_path + 'OSZICAR')
            mag_bare = oszicar.ionic_steps[-1]['mag']
            f.writelines('\tmagnetization: %4.3f\n' % mag_bare)
            df.bare_mag[idx + 1] = float('%4.3f' % mag_bare)

            outcar = Outcar(file_path + 'OUTCAR')
            time_bare = float(outcar.run_stats['Total CPU time used (sec)'])/3600.0
            df.time_bare[idx + 1] = float('%4.3f' % time_bare)

            """
            Magmom list for B sites
            """
            z_param = bare_opt.get_cell_lengths_and_angles()[2]

            metal_idx_list = []
            for atom in bare_opt:
                if atom.symbol != 'O' and atom.position[2]/z_param > 0.29:
                    metal_idx_list.append(atom.index)

            metal_index = metal_idx_list[0]

            metal_magm = []
            for atom in bare_opt:
                if atom.symbol == bare_opt[metal_index].symbol:
                    metal_magm.append(outcar.magnetization[atom.index]['tot'])

            df_mag.magm_bare[idx+1] = str(metal_magm)
                
            write_xsd('results/%03d_%s/bare.xsd' % (idx + 1.0, formula), bare_opt)   # put AFM, def in the path if required
          # view(bare_opt)
 
            """
            Adsorbates
            """
            for k in range(len(adsorbates)):  
                try:
                    if os.path.exists(file_path + '%s/cont/' % adsorbates[k]):
                        file_path_ads = file_path + '%s/cont/' % adsorbates[k]
                    elif os.path.exists(file_path + '%s_np/cont/' % adsorbates[k]):
                        file_path_ads = file_path + '%s_np/cont/' % adsorbates[k]
                    else:
                        file_path_ads = file_path + '%s/' % adsorbates[k]
                    print(file_path_ads)
                   
                    slab_ini = read_vasp(file_path_ads + 'POSCAR')    # put AFM, def in the path if required
                    slab_opt = read_vasp(file_path_ads + 'CONTCAR')   # put AFM, def in the path if required
                    v_slab = Vasprun(file_path_ads + 'vasprun.xml')   # put AFM, def in the path if required

                    E_slab = v_slab.ionic_steps[-1]['e_wo_entrp']
                    DE_ads = E_slab - E_bare - ads_corr[adsorbates[k]]                    
                                          
                    f.writelines('%03d_%s_%s\n' % (idx + 1.0, formula, adsorbates[k]))
                    f.writelines('\tElectronic & ionic converged?: %s\n' % v_slab.converged) 
                    f.writelines('\tFinal_energy(eV): %4.6f\n' % E_slab)
                    f.writelines('\tDE_ads(eV): %4.6f\n' % DE_ads)

                    # new indexing for iloc func.
                    idx_i = idx - num
                    df.iloc[idx_i][7*k + 7] = v_slab.converged
                    df.iloc[idx_i][7*k + 8] = float('%4.4f' % E_slab)
                    df.iloc[idx_i][7*k + 9] = float('%4.3f' % DE_ads)

                    oszicar = Oszicar(file_path_ads + 'OSZICAR')      # put AFM, def in the path if required
                    mag_ads = oszicar.ionic_steps[-1]['mag']
                    dmag = mag_ads - mag_bare
                    
                    f.writelines('\tmagnetization: %4.3f\n' % mag_ads)
                    df.iloc[idx_i][7*k + 10] = float('%4.3f' % mag_ads)
                    df.iloc[idx_i][7*k + 11] = float('%4.3f' % dmag)
                    df.iloc[idx_i][7*k + 12] = len(v_slab.ionic_steps)

                    outcar_ads = Outcar(file_path_ads + 'OUTCAR')      # put AFM, def in the path if required
                    time_ads = float(outcar_ads.run_stats['Total CPU time used (sec)'])/3600.0
                    df.iloc[idx_i][7*k + 13] = float('%4.3f' % time_ads)

                    z_param = slab_opt.get_cell_lengths_and_angles()[2]

                    metal_idx_list = []
                    for atom in slab_opt:
                        if atom.symbol != 'O' and atom.symbol != 'H' and atom.position[2]/z_param > 0.29:
                            metal_idx_list.append(atom.index)

                    if len(metal_idx_list) != 4:
                        f.writelines('\tB_site_indexing error!!\n')

                    metal_index = metal_idx_list[0]
                    metal_magm_ads = []
                    for atom in slab_opt:
                        if atom.symbol == slab_opt[metal_index].symbol:
                            metal_magm_ads.append(outcar_ads.magnetization[atom.index]['tot'])

                    df_mag.iloc[idx_i][k+1] = str(metal_magm_ads)                  
#                   view(slab_ini)
#                   view(slab_opt)

                    print(df)
                    # put AFM, def in the path if required
                    write_xsd('results/%03d_%s/%s_opt.xsd' % (idx + 1.0, formula, adsorbates[k]), slab_opt)        
                  
                except:
                    f.writelines('%03d_%s_%s files are not found !!\n' % (idx + 1.0, formula, adsorbates[k]))                          
        except:
            f.writelines('%03d_%s files are not found !!\n' % (idx + 1.0, formula))
            continue    

    df.to_csv('results/df_summary_%03d_%03d.csv' % (num_ini+1, num_fin))     # put AFM, def in the path if required
    df.to_pickle('results/df_summary_%03d_%03d.pkl' % (num_ini+1, num_fin))  # put AFM, def in the path if required

    df_mag.to_csv('results/df_mag_%03d_%03d.csv' % (num_ini+1, num_fin))     # put AFM, def in the path if required
    df_mag.to_pickle('results/df_mag_%03d_%03d.pkl' % (num_ini+1, num_fin))  # put AFM, def in the path if required

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
