# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 21:27:42 2021

@author: gyjung
"""

import os
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

##############################################################################
def check_convergence(directory):
    try:
        v = Vasprun(directory + 'vasprun.xml')
        if v.converged is True:
            check_results = True
        else:
            check_results = False
    except:
        check_results = 'Unknown'
    return check_results

##############################################################################
# NUPD_dictionary
TM_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
               'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

NUPD_0 = ['Sc','Zn','Y','Cd','Lu','Hg']
NUPD_1 = ['Ti','Cu','Zr','Ag','Hf','Au']
NUPD_2 = ['V','Ni','Nb','Pd','Ta','Pt']
NUPD_3 = ['Cr','Co','Mo','Rh','W','Ir']
NUPD_4 = ['Mn','Fe','Tc','Ru','Re','Os']

NUPD_dict = {}
for TM in TM_elements:
    if TM in NUPD_0:
        NUPD_dict[TM] = [0]
    elif TM in NUPD_1:
        NUPD_dict[TM] = [0,1,2]
    elif TM in NUPD_2:
        NUPD_dict[TM] = [0,1,2,3]
    elif TM in NUPD_3:
        NUPD_dict[TM] = [0,1,2,3,4]
    elif TM in NUPD_4:
        NUPD_dict[TM] = [0,1,2,3,4,5]

properties = ['formula','nupd','convg_fix','tot_E_fix','convg_relax','tot_E_relax','magm_fix','magm_relax','struc_fix','struc_relax']

##############################################################################

initial_time = time.time()

originalPath = os.getcwd()

# Input parameters
model_type = 'dvn4'
# target_elements = ['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd']
target_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
# target_elements = TM_elements


tot_num = 0
for TM in TM_elements:
    tot_num += len(NUPD_dict[TM])

df = pd.DataFrame(np.array([[None for i in range(len(properties))] for j in range(tot_num)]),
                  index = [m for m in range(tot_num)],
                  columns = properties)

num = 0
for idx, TM in enumerate(TM_elements):
    if TM in target_elements:

        if TM == target_elements[0]:
            idx_ini = idx + 1
        if TM == target_elements[-1]:
            idx_fin = idx + 1

        formula = TM + '_' + model_type
        if os.path.exists(originalPath + '/%02d_%s' % (idx+1, formula)):
            os.chdir(originalPath + '/%02d_%s' % (idx+1, formula))
            print(os.getcwd())
        else:
            continue
    
        for nupd in NUPD_dict[TM]:
#           if os.path.exists(os.getcwd() + '/NUPD/opt_%d' % nupd):
            path_1 = os.getcwd() + '/NUPD/opt_%d/' % nupd
        
            for i in range(1,5):
                if os.path.exists(os.getcwd() + '/NUPD/opt_%d/cont%d' % (nupd, i)):
                    path_1 = os.getcwd() + '/NUPD/opt_%d/cont%d/' % (nupd, i)
        
            convg_check_1 = check_convergence(path_1)
            print('%s: %s' % (path_1, convg_check_1))

            path_2 = os.getcwd() + '/NUPD/opt_%d/fin/' % nupd
            convg_check_2 = check_convergence(path_2)
            print('%s: %s' % (path_2, convg_check_2))

#           if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax' % nupd):
            path_3 = os.getcwd() + '/NUPD/opt_%d/relax/' % nupd 

            for j in range(1,5):
                if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax/cont%d' % (nupd, j)):
                    path_3 = os.getcwd() + '/NUPD/opt_%d/relax/cont%d/' % (nupd, j)
        
            convg_check_3 = check_convergence(path_3)
            print('%s: %s' % (path_3, convg_check_3))

            path_4 = os.getcwd() + '/NUPD/opt_%d/relax/fin/' % nupd

            for h in range(1,5):
                if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax/fin/cont%d' % (nupd, h)):
                    path_4 = os.getcwd() + '/NUPD/opt_%d/relax/fin/cont%d/' % (nupd, h)

            convg_check_4 = check_convergence(path_4)
            print('%s: %s' % (path_4, convg_check_4))

            df.formula[num+nupd] = formula
            df.nupd[num+nupd] = nupd

            df.convg_fix[num+nupd] = (convg_check_1, convg_check_2)
            if convg_check_1 is True and convg_check_2 is True:
                calc_1 = Vasp(restart = True, directory = path_2) 
                model_1 = calc_1.get_atoms()           
                df.struc_fix[num + nupd] = model_1
                df.tot_E_fix[num + nupd] = model_1.get_potential_energy()
                df.magm_fix[num + nupd] = model_1.get_magnetic_moment()

                v = Vasprun(path_2 + 'vasprun.xml')
                if v.final_energy > 0:
                    df.tot_E_fix[num+nupd] = v.ionic_steps[-2]['e_wo_entrp']
            
            df.convg_relax[num+nupd] = (convg_check_3, convg_check_4)
            if convg_check_3 is True and convg_check_4 is True:
                calc_2 = Vasp(restart = True, directory = path_4)
                model_2 = calc_2.get_atoms()
                df.struc_relax[num + nupd] = model_2
                df.tot_E_relax[num + nupd] = model_2.get_potential_energy()
                df.magm_relax[num + nupd] = model_2.get_magnetic_moment()

                v = Vasprun(path_4 + 'vasprun.xml')
                E_fin_1 = v.final_energy
                if v.nionic_steps > 1:
                    E_fin_2 = v.ionic_steps[-2]['e_wo_entrp']
                    diff = abs(E_fin_2 - E_fin_1)
                    if diff > 0.1:           
                        df.tot_E_relax[num+nupd] = E_fin_2
            
        num += len(NUPD_dict[TM])

df.to_csv(originalPath + '/results/df_analysis_%02d_%02d.csv' % (idx_ini, idx_fin))
df.to_pickle(originalPath + '/results/df_analysis_%02d_%02d.pkl' % (idx_ini, idx_fin))
               
final_time = time.time()
print('Execution time for script (sec) : %6.1f' % (final_time - initial_time))
 
