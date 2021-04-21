# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:31:14 2021

@author: gyjung
"""

import os
import shutil
import time
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from ase.io import read
from ase.build import surface, make_supercell
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.build import sort        
from ase.calculators.vasp import Vasp

from pymatgen.io.vasp.outputs import Vasprun

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
TM_elements = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
               "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
               "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]

###############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())

###############################################################################
createFolder('models')
df_entries = pd.read_csv('mpid_list_v3.csv')

# insert the number of index for modeling here
num_ini = 309
num_fin = 310

###############################################################################
initial_time = time.time()

for idx in range(num_ini,num_fin):
    start_time = time.time()
    
    formula = df_entries['formula'][idx]

    # Import from bulk results
    if os.path.exists('%03d_%s/2nd/' % (idx + 1.0, formula)):
        if os.path.exists('%03d_%s/2nd/cont/' % (idx + 1.0, formula)):
            file_path = '%03d_%s/2nd/cont/' % (idx + 1.0, formula)
        else:
            file_path = '%03d_%s/2nd/' % (idx + 1.0, formula)
    elif os.path.exists('%03d_%s/2nd_opt/' % (idx + 1.0, formula)):
        if os.path.exists('%03d_%s/2nd_opt/cont/' % (idx + 1.0, formula)):
            file_path = '%03d_%s/2nd_opt/cont/' % (idx + 1.0, formula)
        else:
            file_path = '%03d_%s/2nd_opt/' % (idx + 1.0, formula)

    createFolder(file_path + 'surface')
    
    # Slab modeling
    bulk_opt = read(file_path + 'CONTCAR')
      
    slab = surface(bulk_opt, (1,0,0), 3, vacuum = 7.5)  
    slab_super = make_supercell(slab, [[2,0,0],[0,2,0],[0,0,1]])
    
    positions = slab_super.get_positions()
    
    layer_position = []
    for i in range(len(positions)):
        if positions[i][2] not in layer_position:
            layer_position.append(positions[i][2])
    
    layer_position.sort()

    print('%03d_%s\n    N_atoms(before) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

    natoms_bot = 0
    for atom in slab_super:
        if atom.position[2] < (layer_position[0] + 0.1):
            natoms_bot += 1
        
    if natoms_bot == 8: 
        del slab_super[[atom.position[2] < (layer_position[0] + 0.1) for atom in slab_super]]
        layer_position.remove(layer_position[0])
    elif natoms_bot == 12:
        del slab_super[[atom.position[2] > (layer_position[-1] - 0.1) for atom in slab_super]]
        layer_position.remove(layer_position[-1])
    else:
        print('N_atoms in bottomost layer - %d, which is wrong!!' % natoms_bot)
        print('    N_atoms(after) : %d' % (slab_super.get_global_number_of_atoms()))
    
    fix_id_list = []
    for atom in slab_super:
        if atom.position[2] < (layer_position[-3] + 0.1):
            fix_id_list.append(atom.index)
    
    print('    N_atoms(fixed) : %d' % (len(fix_id_list)))

    slab_super.set_constraint(FixAtoms(indices = fix_id_list))
    
    for atom in slab_super:
        atom.position[2] -= layer_position[0]

    model = sort(slab_super)      
    model.pbc = True
#   view(model)
    elements = model.get_chemical_symbols()    
    
    # Hubbard U setting
    ldau = False
    ldau_luj = {}
    for el in elements:
        if el in U_elements:
            ldau_luj[el] = {'L':2, 'U':U_dict[el], 'J':0}
            ldau = True

    # MAGMOM settings
    mag_dict = {}
    for el in elements:
        if el in TM_elements:
            mag_dict[el] = 4.0
        else:
            mag_dict[el] = 0.6
    magmoms = [mag_dict[el] for el in elements]

    # Spin-polar Geop(1st)
    calc = Vasp(kpts = (4,4,1), system = formula, xc = 'pbe', setups={'base':'materialsproject','W':'_sv'},
                lmaxmix = 4, istart = 0, icharg = 1, encut = 520, ediff = 2e-06, ediffg = -0.02,
                ispin = 2, ibrion = 2, idipol = 3, isif = 2, isym = 0, ismear = 0,
                algo = 'fast', nelm = 500, nelmdl = -12, prec = 'accurate', sigma = 0.05,
                lorbit = 11, npar = 16, lplane = True, ncore = 16, lreal = False, nsw = 200,
                ldau = ldau, ldau_luj = ldau_luj, directory = file_path + 'surface')
    
    model.set_initial_magnetic_moments(magmoms = magmoms)
    model.calc = calc
    try:
        model.get_potential_energy()
        magm = model.get_magnetic_moments()
    except:
         print('    Unknown_error:%s' % file_path)
         continue

    # Spin-polar Geop(2nd)
    createFolder(file_path + 'surface/2nd')
    calc = Vasp(restart = True, directory = file_path + 'surface')
    model_cnt = calc.get_atoms()
    shutil.copy(file_path + 'surface/CHGCAR', file_path + 'surface/2nd/')
    calc.set(ispin = 2, magmom = magm, directory = file_path + 'surface/2nd')
    model_cnt.calc = calc
    try:
        model_cnt.get_potential_energy()
    except:
        print('    Unknown_error:%s' % file_path)

    end_time = time.time()
    print('    Calculation time(sec): %6.1f' % (end_time - start_time))
    
final_time = time.time()
print('Execution time for script (sec) : %6.1f\n' % (final_time - initial_time))
