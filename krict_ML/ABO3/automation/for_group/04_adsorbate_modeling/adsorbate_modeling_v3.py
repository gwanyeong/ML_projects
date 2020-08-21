# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 01:25:53 2020

@author: gyjung
"""

import os
import shutil
import csv
import itertools
import time

from dotenv import load_dotenv

from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp
from ase.visualize import view
from ase.build import add_adsorbate, sort

from pymatgen import MPRester
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
load_dotenv('.env')
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')

mpr = MPRester(MATERIAL_API_KEY)
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

createFolder('models')

# insert the number of index for modeling here
num_ini = 20
num_fin = 30

with open('models/adsorbate_modeling_%03d_%03d.log' % (num_ini + 1, num_fin), 'w') as f:
    start_time = time.time()

    # adsorbate modeling
    OOH = Atoms('OOH',[[0,0,0],[-0.8065, -0.8065, 0.796], [-0.4355, -0.6755, 1.706]])
    O = Atoms('O', [[0,0,0]])
    OH = Atoms('OH', [[0,0,0], [-0.623, -0.623, 0.422]])
    adsorbates = [OOH, O, OH]
    adsorbate_names = ['OOH', 'O', 'OH']

    for idx in range(num_ini, num_fin):
        formula = entries[idx]['pretty_formula']

        try:
            if os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface_np/cont/' % (idx + 1.0, formula)):
                file_path = ('%03d_%s/2nd/surface_np/cont/') % (idx + 1.0, formula)
                f.writelines('(Note) %03d_%s : spin-restricted calc. was initially perforemd\n' % (idx + 1.0, formula))
            elif os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface/2nd/' % (idx + 1.0, formula)):
                file_path = ('%03d_%s/2nd/surface/2nd/') % (idx + 1.0, formula)

            slab = read_vasp(file_path + 'CONTCAR')
            KPOINTS = Kpoints.from_file(file_path + 'KPOINTS')
            potcar_ori = Potcar.from_file(file_path + 'POTCAR')

        except:
            f.writelines('%03d_%s surface result files were not imported\n' % (idx + 1.0, formula))
    
       # view(slab)

        z_param = slab.get_cell_lengths_and_angles()[2]
  
        metal_idx_list = []
    
        for atom in slab:
            if atom.symbol != 'O' and atom.position[2]/z_param > 0.29:
#               print(atom.index," ",atom.symbol," ",atom.position)
                metal_idx_list.append(atom.index)
        
        metal_index = metal_idx_list[0]
        
        f.writelines('%03d_%s_bare\n' % (idx + 1.0, formula))
        f.writelines('\tNumber of %s in top layer: %d\n' % (slab[metal_index].symbol, len(metal_idx_list)))
        f.writelines('\tmetal_index : %d\tposition(x,y) : (%4.3f, %4.3f)\n'
                     % (metal_index, slab[metal_index].position[0],slab[metal_index].position[1]))

        potcar_symbols = potcar_ori.as_dict()['symbols']

        for ads_idx in range(len(adsorbates)):
            INCAR = Incar.from_file(file_path + 'INCAR')
            createFolder(file_path + adsorbate_names[ads_idx])
            add_adsorbate(slab = slab, adsorbate = adsorbates[ads_idx], height = 2.0,
                          position = (slab[metal_index].position[0], slab[metal_index].position[0]))
        
            elements = []
            for element in slab.get_chemical_symbols():
                if element not in elements:
                    elements.append(element)

            for n_atom_ads in range(len(adsorbates[ads_idx])):
                INCAR['MAGMOM'].append(0.6)
        
            mag_dic = dict(zip(slab.get_chemical_symbols(), INCAR['MAGMOM']))
        
            INCAR['MAGMOM'] = [mag_dic[x] for x in sorted(slab.get_chemical_symbols())]  
        
            magm = []
            for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
                magm.append("{}*{}".format(len(tuple(g)), m))
        
            f.writelines('%03d_%s_%s\n' % (idx + 1.0, formula, adsorbate_names[ads_idx]))
            f.writelines('\tformula: %s\n' % (slab.get_chemical_formula()))
            f.writelines(['\tMAGMOM: ', str(magm), '\n'])
                     
        #   for i in range(len(INCAR['MAGMOM'])):
        #       print(sorted(slab.get_chemical_symbols())[i]," ",INCAR['MAGMOM'][i])
               
            slab_sorted = sort(slab)
            del slab[[atom.index > 51 for atom in slab]]
        
            if 'H' in elements:
                INCAR['LDAUL'].append(-1)
                INCAR['LDAUU'].append(0.0)
                INCAR['LDAUJ'].append(0.0)
            
                LDAUL_dic = dict(zip(elements, INCAR['LDAUL']))
                LDAUU_dic = dict(zip(elements, INCAR['LDAUU']))
                LDAUJ_dic = dict(zip(elements, INCAR['LDAUJ']))
            
                INCAR['LDAUL'] = [LDAUL_dic[x] for x in sorted(elements)]
                INCAR['LDAUU'] = [LDAUU_dic[x] for x in sorted(elements)]
                INCAR['LDAUJ'] = [LDAUJ_dic[x] for x in sorted(elements)]
                f.writelines('\tLDAUL: %s\t  LDAUU: %s\tLDAUJ: %s\n' % (INCAR['LDAUL'], INCAR['LDAUU'], INCAR['LDAUJ']))
            
                potcar_symbols.append('H')
                potcar_symbols.sort()
                POTCAR = Potcar(potcar_symbols)
                f.writelines(['\tPOTCAR_symbols: ',str(potcar_symbols), '\n'])
                potcar_symbols.remove('H')
        
            else:
                f.writelines('\tLDAUL: %s\t  LDAUU: %s\tLDAUJ: %s\n' % (INCAR['LDAUL'], INCAR['LDAUU'], INCAR['LDAUJ']))
                POTCAR = Potcar(potcar_symbols)
                f.writelines(['\tPOTCAR_symbols: ',str(potcar_symbols), '\n'])

#           view(slab_sorted)
            write_vasp(file_path + '%s/POSCAR' % (adsorbate_names[ads_idx]), slab_sorted)
           
            INCAR.write_file(file_path + '%s/INCAR' % (adsorbate_names[ads_idx]))
            KPOINTS.write_file(file_path + '%s/KPOINTS' % (adsorbate_names[ads_idx]))     
            POTCAR.write_file(file_path + '%s/POTCAR' % (adsorbate_names[ads_idx]))
        
            # jobscript copy
            destination = file_path + adsorbate_names[ads_idx] + '/'
            job_file = os.getcwd() + '/jobscript_vasp.sh'
            shutil.copy(job_file, destination)  

        f.writelines(['#'*80,'\n'])
    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
