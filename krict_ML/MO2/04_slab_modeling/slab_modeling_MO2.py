# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:40:10 2020

@author: gyjung
"""

import os
import shutil
import itertools
import fileinput
import time
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from ase.io.vasp import read_vasp, write_vasp
from ase.build import surface, make_supercell, sort
from ase.visualize import view
from ase.constraints import FixAtoms

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

##############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

##############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

##############################################################################

TM_elements = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
               "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
               "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]

MO2_list = []
for TM in TM_elements:
    MO2 = TM + "O2"
    MO2_list.append(MO2)


##############################################################################

createFolder('models')

# insert the number of index for modeling here
num_ini = 0
num_fin = len(TM_elements)

##############################################################################

# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
                'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())


# Pseudopotential choice
PP_dict = {'Sc':'Sc_sv', 'Y':'Y_sv', 'Ti':'Ti_pv', 'Zr':'Zr_sv', 'Hf':'Hf_pv',
           'V':'V_sv', 'Nb':'Nb_pv', 'Ta':'Ta_pv', 'Cr':'Cr_pv', 'Mo':'Mo_pv',
           'W':'W_sv', 'Mn':'Mn_pv', 'Tc':'Tc_pv', 'Re':'Re_pv', 'Fe':'Fe_pv',
           'Co':'Co', 'Ni':'Ni_pv', 'Cu':'Cu_pv', 'Zn':'Zn', 'Ru':'Ru_pv',
           'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'Hg':'Hg', 'Au':'Au',
           'Ir':'Ir', 'Pt':'Pt', 'Os':'Os_pv',
           'H':'H', 'O':'O'}

##############################################################################

with open('models/slab_modeling_%03d_%03d.log' % (num_ini + 1, num_fin), 'w') as f:
    start_time = time.time()
    KPOINTS = Kpoints.from_file('KPOINTS_slab')

    for idx in range(num_ini, num_fin):
        formula = MO2_list[idx]
        file_path = '%02d_%s/opt/' % (idx + 1.0, formula)
        createFolder(file_path + 'surface')

        """
        Slab modeling
        """
        bulk_opt = read_vasp(file_path + 'CONTCAR')

        slab = surface(bulk_opt, (1,1,0), 5, vacuum = 7.5)
     #   view(slab)

        slab_super = make_supercell(slab, [[1,0,0],[0,2,0],[0,0,1]])
#       view(slab_super)

        positions = slab_super.get_positions()

        layer_position = []
        for i in range(len(positions)):
            if round(positions[i][2],3) not in layer_position:
                layer_position.append(round(positions[i][2],3))

        layer_position.sort()

        f.writelines('%02d_%s\n\tN_atoms(before) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

        natoms_bot = 0
        for atom in slab_super:
            if atom.position[2] < (layer_position[1] + 0.1):
                natoms_bot += 1

        if natoms_bot == 10:
            del slab_super[[atom.position[2] < (layer_position[1] + 0.1) for atom in slab_super]]
            for layer in layer_position:
                if layer < layer_position[1] + 0.1:
                    layer_position.remove(layer)
        else:
            f.writelines('N_atoms in bottom-most layer - %d, which is wrong!!\n' % natoms_bot)

      #  view(slab_super)

        natoms_top = 0
        for atom in slab_super:
            if atom.position[2] > (layer_position[-1] - 0.1):
                natoms_top += 1

        if natoms_top == 2:
            del slab_super[[atom.position[2] > (layer_position[-1] - 0.1) for atom in slab_super]]
            layer_position.remove(layer_position[-1])
        else:
            f.writelines('N_atoms in topmost layer - %d, which is wrong!!\n' % natoms_top)

     #   view(slab_super)

        f.writelines('    N_atoms(after) : %d' % (slab_super.get_global_number_of_atoms()))

        fix_id_list = []
        for atom in slab_super:
            if atom.position[2] < (layer_position[-7] + 0.1):
                fix_id_list.append(atom.index)

        slab_super.set_constraint(FixAtoms(indices = fix_id_list))
        f.writelines('    N_atoms(fixed) : %d\n' % (len(fix_id_list)))

        slab_super.center()
#       view(slab_super)

        slab_sorted = sort(slab_super)
        view(slab_sorted)


        """
        INCAR
        """
        INCAR = Incar.from_file(file_path + 'INCAR')

        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = 0
        INCAR['IDIPOL'] = 3
        INCAR['NCORE'] = 16

        # Hubbard U Setting
        elements = []
        for element in bulk_opt.get_chemical_symbols():
            if element not in elements:
                elements.append(element)

        LDAUL_dic = dict(zip(elements, INCAR['LDAUL']))
        LDAUU_dic = dict(zip(elements, INCAR['LDAUU']))
        LDAUJ_dic = dict(zip(elements, INCAR['LDAUJ']))

        INCAR['LDAUL'] = [LDAUL_dic[x] for x in sorted(elements)]
        INCAR['LDAUU'] = [LDAUU_dic[x] for x in sorted(elements)]
        INCAR['LDAUJ'] = [LDAUJ_dic[x] for x in sorted(elements)]

        f.writelines(['\t','elements: ',str(sorted(elements)),'\n'])

        """
        MAGMOM
        """
        mag_list = []
        for e in elements:
            if e in TM_elements:
                mag_list.append(5.0)
            else:
                mag_list.append(0.6)

        mag_dic = dict(zip(elements, mag_list))

        magmom = []
        for el in slab_sorted.get_chemical_symbols():
            magmom.append(mag_dic[el])

        INCAR['MAGMOM'] = magmom

        magm = []
        for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
            magm.append("{}*{}".format(len(tuple(g)), m))


        """
        POTCAR
        """
        potcar_ori = Potcar.from_file(file_path + 'POTCAR')
        potcar_symbols = potcar_ori.as_dict()['symbols']
        potcar_symbols.sort()
        POTCAR = Potcar(potcar_symbols)
        f.writelines(['\t','POTCAR_symbols: ', str(potcar_symbols), '\n'])

        f.writelines(['\t', 'LDAUL: ',str(INCAR['LDAUL']),'\t',
                      'LDAUU: ',str(INCAR['LDAUU']),'\t',
                      'LDAUJ: ',str(INCAR['LDAUJ']),'\n'])

        f.writelines(['\t','MAGMOM: ',str(magm),'\n'])

        """
        Write files
        """

        write_vasp(file_path + 'surface/POSCAR', slab_sorted)
        INCAR.write_file(file_path + 'surface/INCAR')
        KPOINTS.write_file(file_path + 'surface/KPOINTS')
        POTCAR.write_file(file_path + 'surface/POTCAR')

        # jobscript copy
        for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
            if '#PBS -N' in line:
                n_line = n
            elif '#PBS -q' in line:
                q_line = n
        PBS_N = '#PBS -N %03d_%s\n' % (idx + 1.0, formula)
        replace_line('jobscript_vasp.sh', n_line, PBS_N)

        queue = 'flat'   # Change queue name if required
        PBS_q = '#PBS -q %s\n' % (queue)
        replace_line('jobscript_vasp.sh', q_line, PBS_q)

        destination = file_path + 'surface/'
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)

        f.writelines(['#'*80,'\n'])

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
     
