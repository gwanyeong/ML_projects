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
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

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

##############################################################################

initial_time = time.time()
originalPath = os.getcwd()

# Input parameters
target_elements = ['Ag', 'Cd']
model_type = 'dvn4'

# path_1 = 'NUPD/opt_%d' or 'NUPD/opt_%d/cont%d'
# path_1_fin = 'NUPD/opt_%d/fin'

# path_2 = 'NUPD/opt_%d/relax' or 'NUPD/opt_%d/relax/cont%d'
# path_2_fin = 'NUPD/opt_%d/relax/fin'

for idx, TM in enumerate(TM_elements):
    if TM in target_elements:

        formula = TM + '_' + model_type
        os.chdir(originalPath + '/%02d_%s' % (idx+1, formula))
        print(os.getcwd())
    
        for nupd in NUPD_dict[TM]:
            if os.path.exists(os.getcwd() + '/NUPD/opt_%d' % nupd):
                path_1 = os.getcwd() + '/NUPD/opt_%d/' % nupd
        
            for i in range(1,5):
                if os.path.exists(os.getcwd() + '/NUPD/opt_%d/cont%d' % (nupd, i)):
                    path_1 = os.getcwd() + '/NUPD/opt_%d/cont%d/' % (nupd, i)
        
            convg_check_1 = check_convergence(path_1)
            print('%s: %s' % (path_1, convg_check_1))

            if convg_check_1 is True:
                start_time = time.time()
                print('    %02d_%s:SPE_nupd_%d_fin' % (idx + 1, formula, nupd))
                createFolder('NUPD/opt_%d/fin' % nupd)
                path_1_fin = os.getcwd() + '/NUPD/opt_%d/fin/' % nupd

                calc_1 = Vasp(restart = True, directory = path_1)
                model_1 = calc_1.get_atoms()
                magm_1 = model_1.get_magnetic_moments()
                model_1.center(vacuum = 10, axis = 2)
                calc_1.set(isif = 2, nelmdl = -12, gamma = True, directory = path_1_fin)
                model_1.set_initial_magnetic_moments(magmoms = magm_1)
                model_1.calc = calc_1

                convg_check_1_fin = check_convergence(path_1_fin)
                if convg_check_1_fin is True:
                    print('%s: %s' % (path_1_fin, convg_check_1_fin))
                else:
                    try:
                        model_1.get_potential_energy()
                    except:
                        print('Unknown error:%02d_%s_fin' % (idx + 1, formula))
                end_time = time.time()
                print('Calculation time(sec) : %6.1f\n' % (end_time - start_time))
              
            if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax' % nupd):
                path_2 = os.getcwd() + '/NUPD/opt_%d/relax/' % nupd 

            for j in range(1,5):
                if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax/cont%d' % (nupd, j)):
                    path_2 = os.getcwd() + '/NUPD/opt_%d/relax/cont%d/' % (nupd, j)
        
            convg_check_2 = check_convergence(path_2)
            print('%s: %s' % (path_2, convg_check_2))

            if convg_check_2 is True:
                start_time = time.time()
                print('    %02d_%s:SPE_nupd_X_fin' % (idx + 1, formula))
                createFolder('NUPD/opt_%d/relax/fin' % nupd)
                path_2_fin = os.getcwd() + '/NUPD/opt_%d/relax/fin/' % nupd

                calc_2 = Vasp(restart = True, directory = path_2)
                model_2 = calc_2.get_atoms()
                magm_2 = model_2.get_magnetic_moments()
                model_2.center(vacuum = 10, axis = 2)
                calc_2.set(isif = 2, nelmdl = -12, gamma = True, directory = path_2_fin)
                model_2.set_initial_magnetic_moments(magmoms = magm_2)
                model_2.calc = calc_2
              
                convg_check_2_fin = check_convergence(path_2_fin)
                if convg_check_2_fin is True:
                    print('%s: %s' % (path_2_fin, convg_check_2_fin))
                else:
                    try:
                        model_2.get_potential_energy()
                    except:
                        print('Unknown error:%02d_%s_fin' % (idx + 1, formula))
                end_time = time.time()
                print('Calculation time(sec) : %6.1f\n' % (end_time - start_time))
               
final_time = time.time()
print('Execution time for script (sec) : %6.1f' % (final_time - initial_time))
 
