# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 11:40:47 2020

@author: gyjung
"""

from ase.io.vasp import read_vasp, write_vasp
from ase.build import surface, make_supercell
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.build import sort        

import os
import shutil

from pymatgen import MPRester
import csv

import itertools

###############################################################################

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating direccotry. ' + directory)


mpr = MPRester('YOUR_MPI_KEY')
      
    
mpid_list = []

###############################################################################
with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243

###############################################################################
entries_from_list = mpr.query(criteria = {"material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula",
                                  "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                  "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                  "e_above_hull","band_gap","volume","theoretical","elements"])
len(entries_from_list)  # 243

entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])


###############################################################################

entries_inc_alkali = mpr.query(criteria = {"elements":{"$all":["O"], "$in":["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                            "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                                            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
                                            "$in":["Sr", "Ba", "K", "Na", "Li", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Si"],    
                                "$nin":["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                                                   "Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
                                                   "Fm","Md","No","Lr"]}, 
                    "material_id":{"$in":mpid_list}},
                    properties = ["material_id","task_id","pretty_formula"])

len(entries_inc_alkali)

entries_alkali_list = []

for i in range(len(entries_inc_alkali)):
    entries_alkali_list.append(entries_inc_alkali[i]['pretty_formula'])

###############################################################################

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar


with open('model_summary_slabs', 'w') as f:
    
     KPOINTS = Kpoints.from_file('KPOINTS_slab')
   
     for idx in range(len(entries)):

        formula = entries[idx]['pretty_formula']
    
        INCAR = Incar.from_file('%03d_%s/2nd/INCAR' % (idx + 1.0, formula))
    
        createFolder('%03d_%s/2nd/surface' % (idx + 1.0, formula))
    
        # INCAR setting
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = 0
        INCAR['IDIPOL'] = 3
        
        """
        Hubbard U setting
        """
        POSCAR = read_vasp('%03d_%s/2nd/CONTCAR' % (idx + 1.0, formula))
    
        elements = POSCAR.get_chemical_symbols()
        for i in range(2):
            elements.remove(elements[-1])
     
        # Correction for previous error
        if INCAR['LDAUL'] == [0.0, 0.0, 0.0]:
            INCAR['LDAUL'] = [-1, -1, -1]
    
        LDAUL_dic = dict(zip(elements, INCAR['LDAUL']))
        LDAUU_dic = dict(zip(elements, INCAR['LDAUU']))
        LDAUJ_dic = dict(zip(elements, INCAR['LDAUJ']))
       
        INCAR['LDAUL'] = [LDAUL_dic[x] for x in sorted(elements)]
        INCAR['LDAUU'] = [LDAUU_dic[x] for x in sorted(elements)]
        INCAR['LDAUJ'] = [LDAUJ_dic[x] for x in sorted(elements)]           
    
    
        # Magmom setting
        if formula in entries_alkali_list:
            INCAR['MAGMOM'] = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                               5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
        else:
            INCAR['MAGMOM'] = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                               5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                               0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
    
        element_list = []
        for n in range(8):
            element_list.append(elements[0])
        for n in range(8,20):
            element_list.append(elements[1])
        for n in range(20,52):
            element_list.append(elements[2])

        mag_dic = dict(zip(element_list, INCAR['MAGMOM']))
        INCAR['MAGMOM'] = [mag_dic[x] for x in sorted(element_list)]
        
        
        magm = []
        for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
            magm.append("{}*{}".format(len(tuple(g)), m))

#        print(sorted(elements))
#        print(INCAR['LDAUL']," ",INCAR['LDAUU']," ",INCAR['LDAUJ'])
#        print(INCAR['MAGMOM'])
        
        f.writelines(['%03d_%s' % (idx + 1.0, formula),'\n'])
        f.writelines(['elements: ',str(sorted(elements)),'\n'])
        
        """
        POTCAR
        """    
        potcar_ori = Potcar.from_file('%03d_%s/2nd/POTCAR' % (idx + 1.0, formula))
        potcar_symbols = potcar_ori.as_dict()['symbols']
        potcar_symbols.sort()
        POTCAR = Potcar(potcar_symbols)
    
#       print(idx + 1.0," ", formula," ",potcar_symbols)
        f.writelines(['POTCAR_symbols: ', str(potcar_symbols), '\n'])

        f.writelines(['LDAUL: ',str(INCAR['LDAUL']),'\t',
                      'LDAUU: ',str(INCAR['LDAUU']),'\t',
                      'LDAUJ: ',str(INCAR['LDAUJ']),'\n'])

        f.writelines(['MAGMOM: ',str(magm),'\n'])

        INCAR.write_file('%03d_%s/2nd/surface/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%03d_%s/2nd/surface/KPOINTS' % (idx + 1.0, formula))     
        POTCAR.write_file('%03d_%s/2nd/surface/POTCAR' % (idx + 1.0, formula))

        # jobscript copy
        destination = '%03d_%s/2nd/surface/' % (idx + 1.0, formula)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)    


        """
        slab modeling
        """    
        bulk_opt = read_vasp('%03d_%s/2nd/CONTCAR' % (idx + 1.0, formula))
      
        slab = surface(bulk_opt, (1,0,0), 3, vacuum = 7.5)
        # view(slab)
    
        slab_super = make_supercell(slab, [[2,0,0],[0,2,0],[0,0,1]])
        #view(slab_super)
    
        positions = slab_super.get_positions()
    
        layer_position = []
        for i in range(len(positions)):
#           print('%4.3f' % positions[i][2],end = '\t')
            if positions[i][2] not in layer_position:
                layer_position.append(positions[i][2])
    
        layer_position.sort()
    
#       view(slab_super)
#       print('%03d_%s - N_atoms(before) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))
        f.writelines('%03d_%s - N_atoms(before) : %d\n' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

        natoms_bot = 0
        for atom in slab_super:
            if atom.position[2] < (layer_position[0] + 0.1):
                natoms_bot += 1
    #   print(natoms_bot)
    
        if natoms_bot == 8: 
            del slab_super[[atom.position[2] < (layer_position[0] + 0.1) for atom in slab_super]]
            layer_position.remove(layer_position[0])
        elif natoms_bot == 12:
            del slab_super[[atom.position[2] > (layer_position[-1] - 0.1) for atom in slab_super]]
            layer_position.remove(layer_position[-1])
        else:
#           print('N_atoms in bottomost layer - %d, which is wrong!! ' % natoms_bot)
            f.writelines('N_atoms in bottomost layer - %d, which is wrong!!\n' % natoms_bot)

        
#        print('%03d_%s - N_atoms(after) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))
        f.writelines('%03d_%s - N_atoms(after) : %d\n' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

    
        fix_id_list = []
        for atom in slab_super:
            if atom.position[2] < (layer_position[-3] + 0.1):
                fix_id_list.append(atom.index)
    
#       print('%03d_%s - N_atoms(fixed) : %d' % (idx + 1.0, formula, len(fix_id_list)))
        f.writelines('%03d_%s - N_atoms(fixed) : %d\n' % (idx + 1.0, formula, len(fix_id_list)))

    
        slab_super.set_constraint(FixAtoms(indices = fix_id_list))
    #   view(slab_super)
    
        for atom in slab_super:
            atom.position[2] -= layer_position[0]

 #      view(slab_super)

        slab_sorted = sort(slab_super)      
        write_vasp('%03d_%s/2nd/surface/POSCAR' % (idx + 1.0, formula), slab_sorted)

        f.writelines(['#'*80,'\n'])

