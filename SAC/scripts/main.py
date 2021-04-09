# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 22:04:24 2021

@author: gyjung
"""

import os
import time
import shutil
import fileinput

from ase import Atom
from ase.io import read, write
from ase.visualize import view
from ase.build import make_supercell, sort
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
        if v.converged:
            check_results = True
        else:
            check_results = False
    except:
        check_results = 'Unknown'
    return check_results

##############################################################################
def recalculation(directory, n_iter):
    current_path = os.getcwd()
    for n, line in enumerate(fileinput.FileInput(directory + 'vasp.out')):
        if 'ZBRENT: fatal error in bracketing' in line:
            print("    %s - restart calc.!" % line)
            os.chdir(directory)
            
            direc_tmp = 'cont%d' % n_iter
            createFolder(direc_tmp)
            calc = Vasp(restart = True)
            atoms = calc.get_atoms()
            calc.set(ispin = 2, directory = direc_tmp)
            atoms.calc = calc
            try:
                energy = atoms.get_potential_energy()
                magmom = atoms.get_magnetic_moments()
            except:
                energy = magmom = None
            os.chdir(current_path)

    return energy, magmom

##############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())
##############################################################################
# NUPD_dictionary
TM_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
               'La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

NUPD_0 = ['Sc','Zn','Y','Cd','La','Hg']
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

target_elements = ['Cu', 'Zn']

# Modeling
for idx, TM in enumerate(TM_elements):
    if TM in target_elements:
    
        formula = TM + '_svn3'
        model = read(originalPath + '/models/svn3_ini.cif')
        model[24].symbol = TM
        model = sort(model)
        write(filename = originalPath + '/models/%02d_%s.cif' % (idx + 1, formula), images = model)

        createFolder(originalPath + '/%02d_%s' % (idx + 1,formula))
        os.chdir(originalPath + '/%02d_%s' % (idx + 1, formula))
        print(os.getcwd())

        elements = model.get_chemical_symbols()
    
        # Hubbard U setting
        ldau = False
        ldau_luj = {}
        for el in elements:
            if el in U_elements:
                ldau_luj[el] = {'L':2, 'U':U_dict[el], 'J':0}
                ldau = True

        calc = Vasp(kpts=(6,5,1), system = formula,
                    xc = 'rpbe', istart = 0, icharg = 1, encut = 520, 
                    ediff = 1e-06, algo = 'fast', # lreal = False
                    ismear = -5, sigma = 0.05, lorbit = 11, 
                    npar = 16, lplane = True, ncore = 16, ldau = ldau, ldau_luj = ldau_luj)    
    
        model.calc = calc
   
        # MAGMOM settings
        mag_dict = {}
        for el in elements:
            if el in TM_elements:
                mag_dict[el] = 4.0
            else:
                mag_dict[el] = 0.6
        magmoms = [mag_dict[el] for el in elements]

        # Spin-restricted SPE
        print('    %02d_%s:spin-restricted SPE.' % (idx + 1, formula))
        createFolder('np')
        calc.set(ispin = 1, nsw = 0, ibrion = -1, directory = 'np')
        try:
            model.get_potential_energy()
        except:
            print('    Unknown error - np:%02d_%s' % (idx + 1, formula))
    
        createFolder('NUPD')
        for nupd in NUPD_dict[TM]:
            start_time = time.time()

            # Spin-polarized Geop
            print('    %02d_%s:geop_nupd_%d' % (idx + 1, formula, nupd))
            createFolder('NUPD/opt_%d' % nupd) 
            model_ini = read('np/POSCAR')
            shutil.copy('np/CHGCAR', 'NUPD/opt_%d/' % nupd)
            calc.set(ispin = 2, nsw = 200, ibrion = 2, isif = 3, ediffg = -0.01, nupdown = nupd,
                     directory = 'NUPD/opt_%d' % nupd)
            model_ini.set_initial_magnetic_moments(magmoms = magmoms)
            model_ini.calc = calc
            try:
                model_ini.get_potential_energy()
                magm = model_ini.get_magnetic_moments()
            except:
                print('    Unknown error:%02d_%s_opt' % (idx + 1, formula))

            # Error check -> restart
            for i in range(1,5):
                if i == 1:
                    direc = 'NUPD/opt_%d/' % nupd
                else:
                    direc = 'NUPD/opt_%d/cont%d/' % (nupd, i-1)

                if check_convergence(direc) == 'Unknown':
                    (energy, magm) = recalculation(direc, i)
                    target_direc = 'NUPD/opt_%d/cont%d/' % (nupd, i)
                elif check_convergence(direc):
                    print('    %s: converged' % direc)
                    target_direc = direc
                    break
            
            # Spin-polarized Geop(2nd)
            print('    %02d_%s:geop_nupd_X' % (idx + 1, formula))
            createFolder('NUPD/opt_%d/relax' % nupd)
            model_cont = read(target_direc + 'CONTCAR')
            shutil.copy(target_direc + 'CHGCAR', 'NUPD/opt_%d/relax/' % nupd)
            calc.set(nupdown = -1, directory = 'NUPD/opt_%d/relax/' % nupd)
            model_cont.set_initial_magnetic_moments(magmoms = magm)
            model_cont.calc = calc
            try:
                model_cont.get_potential_energy()
            except:
                print('    Unknown error:%02d_%s_relax' % (idx + 1, formula))

            # Error check(2nd) -> restart
            for j in range(1,5):
                if j == 1:
                    direc_cnt = 'NUPD/opt_%d/relax/' % nupd
                else:
                    direc_cnt = 'NUPD/opt_%d/relax/cont%d/' % (nupd, j-1)

                if check_convergence(direc_cnt) == 'Unknown':
                    (energy, magm) = recalculation(direc_cnt, j)
                elif check_convergence(direc_cnt):
                    print('    %s: converged' % direc_cnt)
                    break

            end_time = time.time()
            print('      Calc. time(sec): %6.1f' % (end_time - start_time))

final_time = time.time()
print('Execution time for script (sec) : %6.1f' % (final_time - initial_time))
