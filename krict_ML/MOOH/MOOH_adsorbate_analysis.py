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
import warnings
warnings.filterwarnings('ignore')

from dotenv import load_dotenv

from ase.io.vasp import read_vasp, write_vasp
from ase.io.xsd import write_xsd
from ase.visualize import view

from pymatgen.io.vasp import Vasprun, Oszicar, Outcar
from pymatgen import Structure

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################

adsorbates = ['OOH', 'O', 'OH']

H2O_ref = -14.21953641
H2_ref = -6.77110791

OOH_corr = 2*H2O_ref - 1.5*H2_ref
O_corr = H2O_ref - H2_ref
OH_corr = H2O_ref - 0.5*H2_ref

ads_corr = dict(zip(adsorbates,[OOH_corr, O_corr, OH_corr]))

###############################################################################

[OOH_ZPE, OOH_CvdT, OOH_TdS] = [0.46, 0.08, 0.14]
[O_ZPE, O_CvdT, O_TdS] = [0.08, 0.03, 0.05]
[OH_ZPE, OH_CvdT, OH_TdS] = [0.34, 0.05, 0.08]
[H2O_ZPE, H2O_CvdT, H2O_TdS] = [0.57, 0.08, 0.66]
[H2_ZPE, H2_CvdT, H2_TdS] = [0.26, 0.09, 0.40]

G_corr_OOH = OOH_ZPE - (2*H2O_ZPE - 1.5*H2_ZPE) + OOH_CvdT - (2*H2O_CvdT - 1.5*H2_CvdT) -(OOH_TdS - (2*H2O_TdS - 1.5*H2_TdS))
G_corr_O = O_ZPE - (1*H2O_ZPE - H2_ZPE) + O_CvdT -(1*H2O_CvdT - H2_CvdT) -(O_TdS -(1*H2O_TdS - H2_TdS))
G_corr_OH = OH_ZPE - (1*H2O_ZPE - 0.5*H2_ZPE) +OH_CvdT -(1*H2O_CvdT -0.5*H2_CvdT) -(OH_TdS - (1*H2O_TdS - 0.5*H2_TdS))

G_corr_list = [G_corr_OOH, G_corr_O, G_corr_OH]

###############################################################################

MHO2_list = ["ScHO2","TiHO2","VHO2","CrHO2","MnHO2","FeHO2","CoHO2","NiHO2","CuHO2","ZnHO2",
            "YHO2","ZrHO2","NbHO2","MoHO2","TcHO2","RuHO2","RhHO2","PdHO2","AgHO2","CdHO2",
            "HfHO2","TaHO2","WHO2","ReHO2","OsHO2","IrHO2","PtHO2","AuHO2","HgHO2"]
            
createFolder('results/slab101_adsorbate')

properties = ['formula',
              'convg_bare','tot_E_bare','DE_bare','DG_bare','bare_mag','nsteps_bare','struc_bare','time_bare',
              'convg_OOH','tot_E_OOH','DE_OOH','DG_OOH','OOH_mag','OOH_dmag','nsteps_OOH','struc_OOH','time_OOH',
              'convg_O', 'tot_E_O','DE_O', 'DG_O','O_mag','O_dmag','nsteps_O','struc_O','time_O',
              'convg_OH','tot_E_OH','DE_OH','DG_OH','OH_mag','OH_dmag','nsteps_OH','struc_OH','time_OH']

magm_properties = ['magm_bare', 'magm_OOH', 'magm_O', 'magm_OH']

###############################################################################


# insert the number of index for analysis here
num_ini = 0
num_fin = 29
num = num_fin - num_ini

# Determine whether the analysis includes the calc. time and mag from OUTCAR.
read_OUTCAR = True

with open('results/slab101_adsorbate/results_summary_%02d_%02d' % (num_ini+1, num_fin), 'w') as f:     # put AFM, def in the path if required
    start_time = time.time()
   
    df = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                      index = [m for m in range(num_ini+1 ,num_fin+1)],
                      columns = properties)
    df_mag = pd.DataFrame(np.array([['NaN' for p in range(len(magm_properties))] for q in range(num)]),
                          index = [r for r in range(num_ini+1, num_fin+1)],
                          columns = magm_properties)

    for idx, formula in enumerate(MHO2_list):
        df.formula[idx + 1] = formula
        createFolder('results/slab101_adsorbate/%02d_%s' % (idx + 1.0, formula))  # put AFM, def in the path if required
        
        try:
            if os.path.exists(os.getcwd() + '/%02d_%s/cont/slab_101/' % (idx + 1.0, formula)):
                file_path = ('%02d_%s/cont/slab_101/') % (idx + 1.0, formula)
#            elif os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface/2nd/' % (idx + 1.0, formula)):
#                file_path = ('%03d_%s/2nd/surface/2nd/') % (idx + 1.0, formula)
            
            bare_opt = read_vasp(file_path + 'CONTCAR')
            v_bare = Vasprun(file_path + 'vasprun.xml')
            E_bare = v_bare.ionic_steps[-1]['e_wo_entrp']
            struc_bare = Structure.from_file(file_path + 'CONTCAR')

            f.writelines('%02d_%s_bare\n' % (idx + 1.0, formula))
            f.writelines('\tElectronic & ionic converged?: %s\n' % v_bare.converged) 
            f.writelines('\tFinal_energy(eV): %4.6f\n' % E_bare)

            df.convg_bare[idx + 1] = v_bare.converged
            df.tot_E_bare[idx + 1] = float('%4.4f' % E_bare)
            df.DE_bare[idx + 1] = 0
            df.DG_bare[idx + 1] = 0
            df.nsteps_bare[idx + 1] = len(v_bare.ionic_steps)
            df.struc_bare[idx + 1] = bare_opt # struc_bare - Pymatgen Structure object
            
            oszicar = Oszicar(file_path + 'OSZICAR')
            mag_bare = oszicar.ionic_steps[-1]['mag']
            f.writelines('\tmagnetization: %4.3f\n' % mag_bare)
            df.bare_mag[idx + 1] = float('%4.3f' % mag_bare)

            if read_OUTCAR:
                outcar = Outcar(file_path + 'OUTCAR')
                time_bare = float(outcar.run_stats['Total CPU time used (sec)'])/3600.0
                df.time_bare[idx + 1] = float('%4.3f' % time_bare)

                """
                Magmom list for metal
                """
                z_param = bare_opt.get_cell_lengths_and_angles()[2]

                metal_idx_list = []
                for atom in bare_opt:
                    if atom.symbol != 'O' and atom.symbol != 'H':
                        metal_idx_list.append(atom.index)
                print (metal_idx_list)
                metal_index = metal_idx_list[0]

                metal_magm = []
                for atom in bare_opt:
                    if atom.symbol == bare_opt[metal_index].symbol:
                        metal_magm.append(outcar.magnetization[atom.index]['tot'])

                df_mag.magm_bare[idx+1] = str(metal_magm)
            
            write_xsd('results/slab101_adsorbate/%02d_%s/bare.xsd' % (idx + 1.0, formula), bare_opt)   # put AFM, def in the path if required
          # view(bare_opt)
 
            """
            Adsorbates
            """
            for k in range(len(adsorbates)):  
                try:
                    file_path_ad = ('%02d_%s/cont/slab_101/') % (idx + 1.0, formula)
                    if os.path.exists(file_path_ad + '%s/' % adsorbates[k]):
#                        if os.path.exists(file_path + '%s/cont/2nd/' % adsorbates[k]):
#                            file_path_ads = file_path + '%s/cont/2nd/' % adsorbates[k]
#                        else:
#                            file_path_ads = file_path + '%s/cont/' % adsorbates[k]
#                    elif os.path.exists(file_path + '%s_np/cont/' % adsorbates[k]):
#                        if os.path.exists(file_path + '%s_np/cont/2nd/' % adsorbates[k]):
#                            file_path_ads = file_path + '%s_np/cont/2nd/' % adsorbates[k]
#                        else:
#                            file_path_ads = file_path + '%s_np/cont/' % adsorbates[k]
#                    elif os.path.exists(file_path + '%s_np/opt/' % adsorbates[k]):
#                        file_path_ads = file_path + '%s_np/opt/' % adsorbates[k]
#                    else:
                        file_path_ads = file_path_ad + '%s/' % adsorbates[k]
                    print(file_path_ads)
                   
                    slab_ini = read_vasp(file_path_ads + 'POSCAR')    # put AFM, def in the path if required
                    slab_opt = read_vasp(file_path_ads + 'CONTCAR')   # put AFM, def in the path if required
                    v_slab = Vasprun(file_path_ads + 'vasprun.xml')   # put AFM, def in the path if required
                    struc_opt = Structure.from_file(file_path_ads + 'CONTCAR')    # put AFM, def in the path if required

                    E_slab = v_slab.ionic_steps[-1]['e_wo_entrp']
                    DE_ads = E_slab - E_bare - ads_corr[adsorbates[k]]      
                    DG_ads = DE_ads + G_corr_list[k]
                                          
                    f.writelines('%02d_%s_%s\n' % (idx + 1.0, formula, adsorbates[k]))
                    f.writelines('\t%s\n' % file_path_ads)
                    f.writelines('\tElectronic & ionic converged?: %s\n' % v_slab.converged) 
                    f.writelines('\tFinal_energy(eV): %4.6f\n' % E_slab)
                    f.writelines('\tDE_ads(eV): %4.6f\n' % DE_ads)

                    # new indexing for iloc func.
                    idx_i = idx - num_ini
                    df.iloc[idx_i][9*k + 9] = v_slab.converged
                    df.iloc[idx_i][9*k + 10] = float('%4.4f' % E_slab)
                    df.iloc[idx_i][9*k + 11] = float('%4.3f' % DE_ads)
                    df.iloc[idx_i][9*k + 12] = float('%4.3f' % DG_ads)

                    oszicar = Oszicar(file_path_ads + 'OSZICAR')      # put AFM, def in the path if required
                    mag_ads = oszicar.ionic_steps[-1]['mag']
                    dmag = mag_ads - mag_bare
                    
                    f.writelines('\tmagnetization: %4.3f\n' % mag_ads)
                    df.iloc[idx_i][9*k + 13] = float('%4.3f' % mag_ads)
                    df.iloc[idx_i][9*k + 14] = float('%4.3f' % dmag)
                    df.iloc[idx_i][9*k + 15] = len(v_slab.ionic_steps)
                    df.iloc[idx_i][9*k + 16] = slab_opt                # struc_opt (pymatgen Structure object)
 
                    if read_OUTCAR:
                        outcar_ads = Outcar(file_path_ads + 'OUTCAR')      # put AFM, def in the path if required
                        time_ads = float(outcar_ads.run_stats['Total CPU time used (sec)'])/3600.0
                        df.iloc[idx_i][9*k + 17] = float('%4.3f' % time_ads)

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
                    if idx < 10:
                        print(df.head(10))
                    else:
                        print(df[idx-10:idx])
                    # put AFM, def in the path if required
                    write_xsd('results/slab101_adsorbate/%02d_%s/%s_opt.xsd' % (idx + 1.0, formula, adsorbates[k]), slab_opt)        
                  
                except:
                    f.writelines('%02d_%s_%s files are not found !!\n' % (idx + 1.0, formula, adsorbates[k]))
            f.writelines(['#'*29,'\n'])
        except:
            f.writelines('%02d_%s files are not found !!\n' % (idx + 1.0, formula))
            continue    

    df.to_csv('results/slab101_adsorbate/df_summary_%02d_%02d.csv' % (num_ini+1, num_fin))     # put AFM, def in the path if required
    df.to_pickle('results/slab101_adsorbate/df_summary_%02d_%02d.pkl' % (num_ini+1, num_fin))  # put AFM, def in the path if required

    if read_OUTCAR:
        df_mag.to_csv('results/slab101_adsorbate/df_mag_%02d_%02d.csv' % (num_ini+1, num_fin))     # put AFM, def in the path if required
        df_mag.to_pickle('results/slab101_adsorbate/df_mag_%02d_%02d.pkl' % (num_ini+1, num_fin))  # put AFM, def in the path if required

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))