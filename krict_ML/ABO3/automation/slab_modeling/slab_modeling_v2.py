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

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
mpr = MPRester('YOUR_MPI_KEY')
mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

len(mpid_list) # 243

###############################################################################
entries_from_list = mpr.query(criteria = {"material_id":{"$in":mpid_list}},
                    properties = ["pretty_formula","e_above_hull"])
len(entries_from_list)  # 243

entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])


###############################################################################
TM_elements = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
               "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
               "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]

nonmag_elements = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "B", "Si", "O"]

###############################################################################

with open('model_summary_slabs', 'w') as f:
    
     KPOINTS = Kpoints.from_file('KPOINTS_slab')
   
     for idx in range(len(entries)):

        formula = entries[idx]['pretty_formula']   
        createFolder('%03d_%s/2nd/surface' % (idx + 1.0, formula))
    
        """
        Slab modeling
        """    
        bulk_opt = read_vasp('%03d_%s/2nd/CONTCAR' % (idx + 1.0, formula))
      
        slab = surface(bulk_opt, (1,0,0), 3, vacuum = 7.5)  
        slab_super = make_supercell(slab, [[2,0,0],[0,2,0],[0,0,1]])
    
        positions = slab_super.get_positions()
    
        layer_position = []
        for i in range(len(positions)):
#           print('%4.3f' % positions[i][2],end = '\t')
            if positions[i][2] not in layer_position:
                layer_position.append(positions[i][2])
    
        layer_position.sort()
    
        f.writelines('%03d_%s\n\tN_atoms(before) : %d' % (idx + 1.0, formula, slab_super.get_global_number_of_atoms()))

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
            f.writelines('N_atoms in bottomost layer - %d, which is wrong!!\n' % natoms_bot)

        f.writelines('    N_atoms(after) : %d' % (slab_super.get_global_number_of_atoms()))
    
        fix_id_list = []
        for atom in slab_super:
            if atom.position[2] < (layer_position[-3] + 0.1):
                fix_id_list.append(atom.index)
    
        f.writelines('    N_atoms(fixed) : %d\n' % (len(fix_id_list)))

        slab_super.set_constraint(FixAtoms(indices = fix_id_list))
    
        for atom in slab_super:
            atom.position[2] -= layer_position[0]

        slab_sorted = sort(slab_super)      
       #view(slab_sorted)

        """
        INCAR
        """
        INCAR = Incar.from_file('%03d_%s/2nd/INCAR' % (idx + 1.0, formula))
        
        INCAR['ISIF'] = 2
        INCAR['ISMEAR'] = 0
        INCAR['ISYM'] = 0
        INCAR['IDIPOL'] = 3
        
        # Hubbard U setting
        elements = []
        for element in bulk_opt.get_chemical_symbols():
            if element not in elements:
                elements.append(element)
     
        # Correction for previous error
        if INCAR['LDAUL'] == [0.0, 0.0, 0.0]:
            INCAR['LDAUL'] = [-1, -1, -1]
    
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
            if e in nonmag_elements:
                mag_list.append(0.6)
            elif e in TM_elements:
                mag_list.append(5.0)
            else:
                mag_list.append(5.0)

        mag_dic = dict(zip(elements, mag_list))

        magmom = []
        for el in slab_sorted.get_chemical_symbols():
            magmom.append(mag_dic[el])

#       for k in range(len(slab_sorted)):
#           print(slab_sorted[k].symbol," ",magmom[i])
            
        INCAR['MAGMOM'] = magmom

        magm = []
        for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
            magm.append("{}*{}".format(len(tuple(g)), m))

        """
        POTCAR
        """    
        potcar_ori = Potcar.from_file('%03d_%s/2nd/POTCAR' % (idx + 1.0, formula))
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
        write_vasp('%03d_%s/2nd/surface/POSCAR' % (idx + 1.0, formula), slab_sorted)
        INCAR.write_file('%03d_%s/2nd/surface/INCAR' % (idx + 1.0, formula))
        KPOINTS.write_file('%03d_%s/2nd/surface/KPOINTS' % (idx + 1.0, formula))     
        POTCAR.write_file('%03d_%s/2nd/surface/POTCAR' % (idx + 1.0, formula))

        # jobscript copy
        destination = '%03d_%s/2nd/surface/' % (idx + 1.0, formula)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)    

        f.writelines(['#'*80,'\n'])
