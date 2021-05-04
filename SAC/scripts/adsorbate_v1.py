import os
import time
import fileinput
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from ase import Atoms
from ase.io import read
from ase.build import add_adsorbate
from ase.calculators.vasp import Vasp
from ase.visualize import view
from ase.optimize import BFGS

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
#       v = Vasp(restart = True, directory = directory)
        if v.converged is True:
            check_results = True
        else:
            check_results = False
    except:
        check_results = 'Unknown'
    return check_results

##############################################################################
def VaspError(directory):
    check_err = False
    for n, line in enumerate(fileinput.FileInput(directory + 'vasp.out')):
        if 'ZBRENT: fatal error' in line:
            print("    Restart calc. - fatal error")
            check_err = True
        elif "LAPACK: Routine ZPOTRF failed!" in line: # Consider to change 'isym' as 0
#           print("    Restart calc. - LAPACK error")
            check_err = True
        elif "VERY BAD NEWS! internal error in subroutine SGRCON:" in line:
            check_err = 'SGRCON'
    return check_err

##############################################################################
def recalculation(directory, i, isym = 2, copy_chgcar = False): # model
    current_path = os.getcwd() # 'relax/fin/'
    os.chdir(directory)        # 'relax/fin/OOH/'
    direc_tmp = 'cont%d/' % (i)
    createFolder(direc_tmp)

    if check_convergence(direc_tmp) is not True:
        if copy_chgcar is True:
            shutil.copy('CHGCAR', direc_tmp)
        calc_tmp = Vasp()
        calc_tmp.read_json(calc_path + 'settings_init.json')
        calc_tmp.set(isif = 2,directory = direc_tmp, txt = 'vasp.out')
        atoms = read(directory + 'CONTCAR')
        atoms.calc = calc_tmp
    try:
        optimizer = BFGS(atoms = atoms, logfile = 'vasp_opt.log', trajectory = 'images.traj', maxstep = 200)
        optimizer.run(fmax = -0.01)
        energy = atoms.get_potential_energy()
        magmom = atoms.get_magnetic_moments()
    except:
        energy = magmom = None
    os.chdir(current_path)

    return energy, magmom

##############################################################################
def recalculation_loop(directory, n_iter = 4):
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
            if VaspError(target_path) is True:
                (energy, magmom) = recalculation(directory, i) #, model)
            elif VaspError(target_path) == 'SGRCON':
                (energy, magmom) = recalculation(directory, i, isym = -1)
            else:
                unknown_error = True
                break
        elif convg_check is False:
            print('    %s: not converged within given ionic steps' % target_path)
            (energy, magmom) = recalculation(directory, i, copy_chgcar = True)
    if unknown_error is True:
        print("    %s: unknown error - magmom is set to default value" % target_path)
        energy =  None
        magmom = magmoms

    return convg_check, target_path, energy, magmom 

##############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17} # 'W':7.17 for W_sv (WO3) 6.2 for W_pv (WO3)
U_elements = list(U_dict.keys())

##############################################################################
# calc_init = Vasp(xc = 'rpbe', setups = {'base':'materialsproject','W':'_sv'},
#                  kpts = (6,5,1), system = formula, idipol = 3, gamma = True,
#                  istart = 0, icharg = 1, ispin = 2, encut = 520, lmaxmix = 4,
#                  ediff = 1e-06, ediffg = -0.01, algo = 'fast', nelm = 200, nelmdl = -12, 
#                  ibrion = 2, isif = 3, ismear = -5, sigma = 0.05, lorbit = 11, 
#                  npar = 16, lplane = True, ncore = 16,  # isym = 2
#                  txt = 'vasp.out',ldau = ldau, ldau_luj = ldau_luj)

TM_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
               'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

##############################################################################
# input parameters
model_type = 'dvn4'
target_elements = ['Fe']
nupd_calc_dict = {'Sc': 0, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 3, 'Fe': 2, 'Co': 1, 'Ni': 0, 'Cu': 0, 'Zn': 0}

# adsorbate modeling
OOH = Atoms('OOH',[[0,0,0],[-1.372,-0.156,0.441],[-1.277,-0.510,1.363]])
O = Atoms('O',[[0,0,0]])
OH = Atoms('OH',[[0,0,0],[-0.858,-0.025,0.562]])
Cl = Atoms('Cl',[[0,0,0]])

adsorbates = [OOH,O,OH,Cl]
adsorbate_names = ['OOH','O','OH','Cl']

##############################################################################
initial_time = time.time()
originalPath = os.getcwd()

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
        calc_path = originalPath + '/%02d_%s/' % (idx + 1, formula)

        nupd = nupd_calc_dict[TM]
        slab_path = originalPath + '/%02d_%s/NUPD/opt_%d/relax/fin/' % (idx+1, formula, nupd)
        for i in range(1,5):
            if os.path.exists(slab_path + 'cont%d' % i):
                slab_path += 'cont%d/' % i

        if check_convergence(slab_path) is not True:
            print("%s: not converged" % slab_path)
            continue

        slab = read(slab_path + 'CONTCAR')
#       print(slab[-1].symbol)
        n_el = len(slab)

        os.chdir(slab_path)

        for ads_idx in range(len(adsorbates)):
            start_time = time.time()
            add_adsorbate(slab = slab, adsorbate = adsorbates[ads_idx], height = 2.0,
                          position = (slab[-1].position[0], slab[-1].position[1]))
#           view(slab)
            adsorbate = adsorbate_names[ads_idx]
            elements = slab.get_chemical_symbols()

            # Hubbard U setting
            ldau = False
            ldau_luj = {}
            for el in elements:
                if el in U_elements:
                    ldau_luj[el] = {'L':2,'U':U_dict[el],'J':0}
                    ldau = True

            # Magmom setting
            calc_bare = Vasp(restart = True, directory = slab_path)
            magmoms = calc_bare.get_magnetic_moments()
            for atom in adsorbates[ads_idx]:
                magmoms = np.append(magmoms, 0.6)
            print(magmoms)

            # Spin-polar geop (1st)
            createFolder(adsorbate)
            if check_convergence(adsorbate + '/') is not True:
                calc_1 = Vasp()
                calc_1.read_json(calc_path + 'settings_init.json')
                calc_1.set(ldau = ldau, ldau_luj = ldau_luj, isif = 2, magmom = magmoms,
                           directory = adsorbate, txt = 'vasp.out')
                slab.calc = calc_1
                optimizer = BFGS(atoms = slab, logfile = adsorbate + '/vasp_opt.log', maxstep = 200,
                                 trajectory = adsorbate + '/images.traj')
                try:
                    optimizer.run(fmax = -0.02) 
                except:
                    print('    Unknown error:%02d_%s_%s' % (idx + 1, formula, adsorbate))

            # Error check -> restart
            (convg_check_1, path_1, energy_1, magm_1) = recalculation_loop(directory = adsorbate + '/') 

            # Spin-polar geop (2nd)
            createFolder(adsorbate + '/2nd')
            if check_convergence(adsorbate + '/2nd/') is not True:
                slab_cnt = read(path_1 + 'CONTCAR')
                calc_2 = Vasp()
                calc_2.read_json(calc_path + 'settings_init.json')
                calc_2.set(ldau = ldau, ldau_luj = ldau_luj, isif = 2, magmom = magm_1,
                           directory = adsorbate + '/2nd', txt = 'vasp.out')
                slab_cnt.calc = calc_2
                optimizer = BFGS(atoms = slab_cnt, logfile = adsorbate + '/2nd/vasp_opt.log', maxstep = 200,
                                 trajectory = adsorbate + '/2nd/images.traj') 
                try:
                    optimizer.run(fmax = -0.02)
                except:
                    print('    Unknown error:%02d_%s_%s' % (idx + 1, formula, adsorbate))

            # Error check -> restart
            (convg_check_2, path_2, energy_2, magm_2) = recalculation_loop(directory = adsorbate + '/2nd/')

            del slab[[atom.index >= n_el for atom in slab]]
            end_time = time.time()
            print('    Calculation time(sec) : %6.1f' % (end_time - start_time))

final_time = time.time()
print('Execution time for script (sec) : %6.1f' % (final_time - initial_time))

