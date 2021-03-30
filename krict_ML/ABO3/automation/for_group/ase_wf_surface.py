# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:45:33 2021

@author: gyjung
"""

import os
import shutil
import time
import pandas as pd

from ase.io.vasp import read_vasp
from ase.build import surface, make_supercell, sort
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp

# from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.io.vasp import Vasprun

##############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)
#############################################################################
def check_convergence(directory):
    v = Vasprun(directory + '/vasprun.xml')
    if not v.converged:
        print('Error: not converging!!')
          
##############################################################################
TM_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
               'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())

##############################################################################

# createFolder('models')

df_entries = pd.read_csv('formula_list.csv')

# insert the number of index for modeling here
num_ini = 0
num_fin = 10

##############################################################################
initial_time = time.time()

for idx in range(num_ini, num_fin):
    start_time = time.time()

    formula = df_entries['formula'][idx]
    file_path = '%03d_%s/' % (idx + 1, formula)
    createFolder(file_path + 'surface')
    slab_path = file_path + 'surface/'

    # Slab modeling
    try:
        bulk_opt = read_vasp(file_path + '2nd/CONTCAR')
    except:
        print("File not found!!: %d_%s\n" % (idx + 1, formula))
        continue

    slab = surface(bulk_opt, (1,0,0), 3, vacuum = 7.5, periodic = True)
    slab_super = make_supercell(slab, [[2,0,0],[0,2,0],[0,0,1]])
    
    positions = slab_super.get_positions()
    
    layer_position = []
    for i in range(len(positions)):
        if positions[i][2] not in layer_position:
            layer_position.append(positions[i][2])
            
    layer_position.sort()
    
    print('%03d_%s\n\tN_atoms(before) : %d' \
          % (idx + 1, formula, slab_super.get_global_number_of_atoms()))
    
    natoms_bot = 0
    for atom in slab_super:
        if atom.position[2] < (layer_position[0] + 0.1):
            natoms_bot += 1
    if natoms_bot == 8:
        del slab_super[[atom.position[2] < (layer_position[0]+0.1) for atom in slab_super]]
        layer_position.remove(layer_position[0])
    elif natoms_bot == 12:
        del slab_super[[atom.position[2] > (layer_position[-1]-0.1) for atom in slab_super]]
        layer_position.remove(layer_position[-1])
    else:
        print('N_atoms in bottommost layer - %d, which is wrong!!\n' % natoms_bot)
    
    print('    N_atoms(after) : %d' % slab_super.get_global_number_of_atoms())
    
    fix_id_list = []
    for atom in slab_super:
        if atom.position[2] < (layer_position[-3] + 0.1):
            fix_id_list.append(atom.index)
    
    print('    N_atoms(fixed) : %d\n' % len(fix_id_list))
    
    slab_super.set_constraint(FixAtoms(indices = fix_id_list))
    
    for atom in slab_super:
        atom.position[2] -= layer_position[0]
        
    model = sort(slab_super)
#   view(model)
    elements = model.get_chemical_symbols()

    # MAGMOM settings
    mag_dict ={}
    for el in elements:
        if el in TM_elements:
            mag_dict[el] = 5.0  # ferromagnetic
        else:
            mag_dict[el] = 0.6
#   print(mag_dict)
    magmoms = [mag_dict[el] for el in elements]      # Edit here if AFM or NM

    # Hubbard U setting   
    ldau = False
    ldau_luj = {}
    for el in elements:
        if el in U_elements:
            ldau_luj[el] = {'L':2,'U':U_dict[el],'J':0}
            ldau = True

    # Will be updated to use 'json' file instead.
    calc = Vasp(kpts = (4,4,1), system = formula + '_slab',
                xc = 'pbe', istart = 0, icharg = 1, encut = 520, 
                ediff = 2e-06, lreal = False, algo = 'fast', ediffg = -0.02,
                nelmdl = -12, nelm = 500, prec = 'accurate', isif = 2,
                ismear = 0, sigma = 0.05, idipol = 3, 
                npar = 16, lplane = True, ncore = 16, ldau = ldau, ldau_luj = ldau_luj)
    
    model.calc = calc

    # Spin-restricted SPE
    print('    Spin-restricted SPE calc.')
    createFolder(slab_path + 'np')
    calc.set(ispin = 1, nsw = 0, ibrion = -1, directory = slab_path + 'np')
    try:
        model.get_potential_energy()
    except:
        print('Error!')
        continue
        
    # spin-unrestricted Geop
    print('    Spin-polar geop.:')
    createFolder(slab_path + 'opt')
    shutil.copy(slab_path + 'np/CHGCAR', slab_path + 'opt/')
    calc.set(ispin = 2, nsw = 200, ibrion = 2, directory = slab_path + 'opt')
    model.set_initial_magnetic_moments(magmoms = magmoms)
    try:
        model.get_potential_energy()
        check_convergence('opt')
    except:
        print('Error!')
        continue
    
    # spin-unrestricted Geop (2nd)
    print('    Spin-polar geop.(2nd)')
    createFolder(slab_path + '2nd')
    shutil.copy(slab_path + 'opt/CHGCAR', slab_path + '2nd/')
    calc.set(ispin = 2, nsw = 200, ibrion = 2, directory = slab_path + '2nd')
    try:
        model.get_potential_energy()
        check_convergence('2nd')
    except:
        print('Error!')
        continue

    end_time = time.time()
    print('    Calc. time(sec): %6.1f' % (end_time - start_time))
    
final_time = time.time()
print('Total execution time for script (sec): %6.1f\n' % (final_time - initial_time))    

