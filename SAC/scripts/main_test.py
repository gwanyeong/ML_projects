# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 21:27:42 2021

@author: gyjung
"""

import os
import time
import shutil
import fileinput
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from ase.io import read, write
# from ase.visualize import view
from ase.build import sort
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
def set_vacuum(model, TM):
    z_length = model.get_cell_lengths_and_angles()[2]
    for atom in model:
        if atom.position[2]/z_length > 0.5:
            atom.position[2] -= z_length
        if atom.symbol == TM and atom.position[2] < 0:
            atom.position[2] = -1 * atom.position[2]
    model.center(vacuum = 10, axis = 2)

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
def VaspFatalError(directory):
    check_err = False
    for n, line in enumerate(fileinput.FileInput(directory + 'vasp.out')):
        if 'ZBRENT: fatal error in bracketing' in line:
            print("    Restart calc. - %s" % line)
            check_err = True
    return check_err

##############################################################################
def recalculation(directory, i, fname, copy_chgcar = False):
    current_path = os.getcwd() # 'NUPD/opt_%d/' % nupd
    os.chdir(directory)
    direc_tmp = 'cont%d' % (i)
    createFolder(direc_tmp)

    if copy_chgcar is True:
        shutil.copy('CHGCAR', direc_tmp + '/')
    calc_tmp = Vasp()
    calc_tmp.read_json(directory + fname)
    calc_tmp.set(ispin = 2, directory = direc_tmp, txt = 'vasp.out')
    atoms = read(directory + 'CONTCAR')
    atoms.calc = calc_tmp
    try:
        energy = atoms.get_potential_energy()
        magmom = atoms.get_magnetic_moments()
    except:
        energy = magmom = None
    os.chdir(current_path)

    return energy, magmom

##############################################################################
def recalculation_loop(directory, n_iter, fname):
    unknown_error = None
    for i in range(1, n_iter + 1):
        if i == 1:
            target_path = directory # 'NUPD/opt_%d/' % nupd
        else:
            target_path = directory + 'cont%d/' % (i-1)

        convg_check = check_convergence(target_path)
        if convg_check is True:
            print('    %s: converged' % target_path)
            calc_tmp = Vasp(restart = True, directory = target_path)
            atoms = calc_tmp.get_atoms()
            energy = atoms.get_potential_energy()
            magmom = atoms.get_magnetic_moments()
            break
        elif convg_check == 'Unknown':
            if VaspFatalError(target_path) is True:
                (energy, magm) = recalculation(directory, i, fname)
            else:
                unknown_error = True
                break
        elif convg_check is False:
            print('    %s: not converged within given ionic steps' % target_path)
            (energy, magm) = recalculation(directory, i, fname, copy_chgcar = True)
    if unknown_error is True:
        energy = magmom = None
        print('    %s: unknown error!' % target_path)
    return convg_check, target_path, energy, magmom 

##############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17} # 'W':7.17 for W_sv (WO3) 6.2 for W_pv (WO3)
U_elements = list(U_dict.keys())

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
target_elements = ['Cu', 'Zn']
model_type = 'svc3'

# path_1 = 'NUPD/opt_%d' 
# path_1_cnt = 'NUPD/opt_%d' or 'NUPD/opt_%d/cont%d'
# path_1_fin = 'NUPD/opt_%d/fin'

# path_2 = 'NUPD/opt_%d/relax' 
# path_2_cnt = 'NUPD/opt_%d/relax' or 'NUPD/opt_%d/relax/cont%d'
# path_2_fin = 'NUPD/opt_%d/relax/fin'

for idx, TM in enumerate(TM_elements):
    if TM in target_elements:

        formula = TM + '_' + model_type

        # Modeling
        model = read(originalPath + '/models/%s_ini.cif' % model_type)
        if model_type.startswith('dv'):
            model[-1].symbol = TM
        elif model_type.startswith('sv'):
            model[24].symbol = TM
        model = sort(model)
        write(filename = originalPath + '/models/%02d_%s.cif' % (idx + 1, formula), images = model)

        createFolder(originalPath + '/%02d_%s' % (idx + 1, formula))
        elements = model.get_chemical_symbols()

        # Hubbard U setting
        ldau = False
        ldau_luj = {}
        for el in elements:
            if el in U_elements:
                ldau_luj[el] = {'L':2,'U':U_dict[el],'J':0}
                ldau = True

        # MAGMOM settings
        mag_dict = {}
        for el in elements:
            if el in TM_elements:
                mag_dict[el] = 4.0
            else:
                mag_dict[el] = 0.6
        magmoms = [mag_dict[el] for el in elements]

        os.chdir(originalPath + '/%02d_%s' % (idx+1, formula))
        print(os.getcwd())

        # Calc. setting
        calc_init = Vasp(xc = 'rpbe', setups = {'base':'materialsproject','W':'_sv'},
                         kpts = (6,5,1), system = formula, idipol = 3, gamma = True,
                         istart = 0, icharg = 1, encut = 520, lmaxmix = 4,
                         ediff = 1e-06, algo = 'fast', nelm = 200, nelmdl = -12,
                         ismear = -5, sigma = 0.05, lorbit = 11,
                         npar = 16, lplane = True, ncore = 16, ldau = ldau, ldau_luj = ldau_luj)

        calc_init.write_json('settings_init.json')
        
        if os.path.exists('np'):
            print('    %02d_%s:spin-restricted SPE: Done' % (idx + 1, formula))

        else:
            # Spin-restricted SPE
            print('    %02d_%s:spin-restricted SPE' % (idx + 1, formula))
            createFolder('np')
            calc = Vasp()
            calc.read_json('settings_init.json')
            calc.set(ispin = 1, nsw = 0, ibrion = -1, directory = 'np', txt = 'vasp.out')
            model.calc = calc
            try:
                model.get_potential_energy()
            except:
                print('    Unknown error - np:%02d_%s' % (idx + 1, formula))

        createFolder('NUPD')
        for nupd in NUPD_dict[TM]:
            start_time = time.time()
            filename_1 = 'settings_nupd_%d.json' % nupd
            if os.path.exists(os.getcwd() + '/NUPD/opt_%d' % nupd):
                path_1 = os.getcwd() + '/NUPD/opt_%d/' % nupd
            else:
                createFolder('NUPD/opt_%d' % nupd)
                model_1 = read('np/POSCAR')
                elements = model_1.get_chemical_symbols()

                shutil.copy('np/CHGCAR', 'NUPD/opt_%d/' % nupd)
                calc_1 = Vasp()
                calc_1.read_json('settings_init.json')
                calc_1.set(ispin = 2, nsw = 200, ibrion = 2, isif = 3, ediffg = -0.01, nupdown = nupd,
                           directory = 'NUPD/opt_%d' % nupd, txt = 'vasp.out')
                model_1.set_initial_magnetic_moments(magmoms = magmoms)
                model_1.calc = calc_1
                calc_1.write_json(filename_1)

                try:
                    model_1.get_potential_energy()
                    magm_1 = model_1.get_magnetic_moments()
                except:
                    print('    Unknown error:%02d_%s' % (idx + 1, formula))
                path_1 = os.getcwd() + '/NUPD/opt_%d/' % nupd

            # Error check -> restart
            (convg_check_1, path_1_cnt, energy_1, magm_1) = recalculation_loop(directory = path_1, n_iter = 4, fname = filename_1)
            print('    %s: %s' % (path_1_cnt, convg_check_1))

            if convg_check_1 is True:
                print('    %02d_%s_nupd_%d_fin' % (idx + 1, formula, nupd))
                createFolder('NUPD/opt_%d/fin' % nupd)
                path_1_fin = os.getcwd() + '/NUPD/opt_%d/fin/' % nupd

                calc_1 = Vasp(restart = True, directory = path_1_cnt)
                model_1 = calc_1.get_atoms()
                magm_1 = model_1.get_magnetic_moments()

#                calc_1 = Vasp()
#                calc_1.read_json('NUPD/opt_%d/settings_nupd_%d.json' % (nupd,nupd))
#                model_1 = read(path_1_cnt + 'CONTCAR')
                
                set_vacuum(model_1, TM)
                calc_1.set(isif = 2, magmom = magm_1, nelmdl = -12, gamma = True, isym = 2, directory = path_1_fin, txt = 'vasp.out')
                model_1.calc = calc_1
                filename_1_fin = 'settings_nupd_%d_fin.json' % nupd
                calc_1.write_json(filename_1_fin)

                convg_check_1_fin = check_convergence(path_1_fin)
                if convg_check_1_fin is True:
                    print('    %s: %s' % (path_1_fin, convg_check_1_fin))
                else:
                    try:
                        model_1.get_potential_energy()
                    except:
                        print('Unknown error:%02d_%s_fin' % (idx + 1, formula))
                    (convg_check_1_fin, path_1_fin, energy_1, magm_1) = recalculation_loop(directory = path_1_fin, n_iter = 4, fname = filename_1_fin)
                    print('    %s: %s' % (path_1_fin, convg_check_1_fin))

            end_time = time.time()
            print('Calculation time(sec) : %6.1f\n' % (end_time - start_time))

            # Spin-polar geop(NUPD_free)
            start_time = time.time()
            filename_2 = 'settings_nupd_X.json'
            if os.path.exists(os.getcwd() + '/NUPD/opt_%d/relax' % nupd):
                path_2 = os.getcwd() + '/NUPD/opt_%d/relax/' % nupd 
            else:
                createFolder('NUPD/opt_%d/relax' % nupd)
                path_2 = os.getcwd() + '/NUPD/opt_%d/relax/' % nupd
                if convg_check_1 is True:
                    shutil.copy(path_1_cnt + 'CHGCAR', path_2)
                model_2 = read(path_1_cnt + 'CONTCAR')

                calc_2 = Vasp()
                calc_2.read_json(path_1 + filename_1)

                calc_2.set(magmom = magm_1, nelmdl = -12, nupdown = -1, gamma = True, isym = 2, directory = path_2, txt = 'vasp.out')
                model_2.calc = calc_2
                calc_2.write_json(filename_2)
                try:
                    model_2.get_potential_energy()
                except:
                    print('    Unknown error! - NUPD/opt_%d/relax' % nupd)

            # Error check -> restart
            (convg_check_2, path_2_cnt, energy_2, magm_2) = recalculation_loop(directory = path_2, n_iter = 4, fname = filename_2)
            print('    %s: %s' % (path_2_cnt, convg_check_2))

            if convg_check_2 is True:
                print('    %02d_%s_nupd_X_fin' % (idx + 1, formula))
                createFolder('NUPD/opt_%d/relax/fin' % nupd)
                path_2_fin = os.getcwd() + '/NUPD/opt_%d/relax/fin/' % nupd
               
                calc_2 = Vasp(restart = True, directory = path_2_cnt)
                model_2 = calc_2.get_atoms()
                magm_2 = model_2.get_magnetic_moments()

#               calc_2 = Vasp()
#               calc_2.read_json('NUPD/opt_%d/relax/settings_nupd_X.json' % nupd)
#               model_2 = read(path_2_cnt + 'CONTCAR')

                set_vacuum(model_2, TM)
                calc_2.set(isif = 2, magmom = magm_2, nelmdl = -12, gamma = True, isym = 2, directory = path_2_fin, txt = 'vasp.out')
                model_2.calc = calc_2
                filename_2_fin = 'settings_nupd_X_fin.json'
                calc_2.write_json(filename_2_fin)
              
                convg_check_2_fin = check_convergence(path_2_fin)
                if convg_check_2_fin is True:
                    print('    %s: %s' % (path_2_fin, convg_check_2_fin))
                else:
                    try:
                        model_2.get_potential_energy()
                    except:
                        print('Unknown error:%02d_%s_fin' % (idx + 1, formula))
                    (convg_check_2_fin, path_2_fin, energy_2, magm_2) = recalculation_loop(directory = path_2_fin, n_iter = 4, fname = filename_2_fin)
                    print('    %s: %s' % (path_2_fin, convg_check_2_fin))

            end_time = time.time()
            print('Calculation time(sec) : %6.1f\n' % (end_time - start_time))

final_time = time.time()
print('Execution time for script (sec) : %6.1f' % (final_time - initial_time))
 
