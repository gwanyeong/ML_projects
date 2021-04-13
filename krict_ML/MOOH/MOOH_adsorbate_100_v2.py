# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 01:25:53 2020
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

from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp
from ase.visualize import view
from ase.build import add_adsorbate, sort

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Vasprun

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

###############################################################################
MOOH_list = ["ScHO2","TiHO2","VHO2","CrHO2","MnHO2","FeHO2","CoHO2","NiHO2","CuHO2","ZnHO2",
            "YHO2","ZrHO2","NbHO2","MoHO2","TcHO2","RuHO2","RhHO2","PdHO2","AgHO2","CdHO2",
            "HfHO2","TaHO2","WHO2","ReHO2","OsHO2","IrHO2","PtHO2","AuHO2","HgHO2"]

#v4 modeling target list
#cont_list = ["VH2O2","CoH2O2","NiH2O2","CdH2O2","OsH2O2"]
#const_list = ["ZrH2O2","MoH2O2","TcH2O2","RhH2O2","PdH2O2","AgH2O2","HfH2O2","WH2O2","PtH2O2","AuH2O2"]
###############################################################################
createFolder('models_ads')
with open('models_ads/adsorbate_modeling.log','w') as f:
    start_time=time.time()

    #adsorbate modeling
    OOH = Atoms('OOH',[[0,0,0],[0,-0.9159,-1.0804],[0,-0.319,-1.8594]]) #left-high / right-low
    OOH2 = Atoms('OOH',[[0,0,0],[0,0.9159,-1.0804],[0,0.319,-1.8594]]) #left-low / right-high
    O = Atoms('O',[[0,0,0]])
    OH = Atoms('OH',[[0,0,0],[0,-0.98,-0.0]])
    OH2 = Atoms('OH',[[0,0,0],[0,0.5296,-0.8171]])
    adsorbates = [OOH, O, OH]
    adsorbates = [OOH, O, OH]
    adsorbate_names = ['OOH','O','OH']

    for idx, formula in enumerate(MOOH_list):
        try:
        
            f.writelines('%02d_%s\n' % (idx + 1.0,formula)) 
            TM = formula.replace("HO2","")
         
            file_path = ('%02d_%s/cont/slab_100/cont5/cont6/cont7/') % (idx + 1.0,formula)
            new_path = ('%02d_%s/cont/slab_100/') % (idx + 1.0,formula)
            createFolder (new_path + 'adsorbate')
            
            slab = read_vasp(file_path + 'CONTCAR')
            KPOINTS = Kpoints.from_file(file_path + 'KPOINTS')
            potcar_ori = Potcar.from_file(file_path + 'POTCAR')
            potcar_symbols = potcar_ori.as_dict()['symbols']
            v = Vasprun(file_path + 'vasprun.xml')

            f.writelines('(File_path) %s\n' % file_path)
            f.writelines('%02d_%s (100) surface - Electronic & ionic converged? %s\n' % (idx + 1.0, formula, v.converged))
            
            metal_z_list = []
            O_z_list = []
            atom_z_list = []
            
            for atom in slab:
                atom_z_list.append(atom.position[2])
                if atom.symbol != 'O' and atom.symbol != 'H':
                    metal_z_list.append(atom.position[2])
                elif atom.symbol == 'O':
                    O_z_list.append(atom.position[2])            
            metal_z_list.sort()
            O_z_list.sort()
            atom_z_list.sort()
            O_z_list
            
            for atom in slab:
                if atom.symbol != 'O' and atom.symbol != 'H' and atom.position[2] == metal_z_list[0]:
                    metal_index = atom.index
                    
                elif atom.symbol != 'O' and atom.symbol != 'H' and atom.position[2] == metal_z_list[2]: 
                    metal_2nd_index = atom.index
                
                elif atom.symbol == 'O' and atom.position[2] == O_z_list[0]:
                    O_index = atom.index
                    
                elif atom.symbol == 'O' and atom.position[2] == O_z_list[-1]:
                    O_index2 = atom.index
                    
                elif atom.symbol == 'O' and atom.position[2] == O_z_list[3]:
                    O_2nd_index = atom.index
                    
                elif atom.symbol == 'O' and atom.position[2] == O_z_list[5]:
                    O_3rd_index = atom.index
                    
                elif atom.position[2] == atom_z_list[0]:
                    top_index = atom.index
                    
            
            #v4 x-axis exception case correction
#            if slab.get_distance(metal_index,O_index) < slab.get_distance(metal_index,O_index2):
#                O_index_real = O_index
#                if abs(slab[metal_index].position[0] - slab[O_index].position[0]) > 1:
#                    O_index_real = O_index2
#            elif slab.get_distance(metal_index,O_index) > slab.get_distance(metal_index,O_index2):
#                O_index_real = O_index2
#                if abs(slab[metal_index].position[0] - slab[O_index2].position[0]) > 1:
#                    O_index_real = O_index
            
            f.writelines('%02d_%s_slab_optimized\n' % (idx + 1.0, formula))
            f.writelines('\tTop metal position(x,y,z): (%4.3f, %4.3f, %4.3f)\n'
                         % (slab[metal_index].position[0], slab[metal_index].position[1], slab[metal_index].position[2]))
            f.writelines('\tSymmetric O position(x,y,z): (%4.3f, %4.3f, %4.3f)\n'
                         % (slab[O_index].position[0], slab[O_index].position[1], slab[O_index].position[2]))
            f.writelines('\tTop atom position(%s)(x,y,z): (%4.3f, %4.3f, %4.3f)\n' % (slab[top_index].symbol,
                            slab[top_index].position[0], slab[top_index].position[1], slab[top_index].position[2]))

            cell_y = slab.get_cell_lengths_and_angles()[1]
            
            disp_z = abs(slab[metal_2nd_index].position[2] - slab[O_2nd_index].position[2])
            disp_y = - abs(slab[metal_2nd_index].position[1] - slab[O_2nd_index].position[1])
            disp_y2 = - abs(slab[metal_2nd_index].position[1] - slab[O_2nd_index].position[1] + cell_y)
            disp_z2 = abs(slab[metal_index].position[2] - slab[O_3rd_index].position[2])
            disp_y3 = - abs(slab[metal_index].position[1] - slab[metal_2nd_index].position[1])
            
            z_correc = slab[metal_index].position[2] - slab[O_index2].position[2]
            print (disp_z)
#            dist = slab.get_distance(metal_index,O_index)
#            z_coord = slab[metal_index].position[2] + disp_z
            
#            if disp_y/cell_y < 0.5:
#                disp_y = cell_y - disp_y
            
            #adsorbate+slab modeling
            for ads_idx in range(len(adsorbates)):
                try:
                    INCAR = Incar.from_file(file_path + 'INCAR')
                    file_write_path = ('%02d_%s/cont/slab_100/adsorbate/') % (idx + 1.0,formula)
                    createFolder(file_write_path + '/' + adsorbate_names[ads_idx])
                    
                    if disp_z < 2.5:
                                        
                        if slab[O_2nd_index].position[1]/cell_y < 0.5:
                            add_adsorbate(slab = slab, adsorbate = adsorbates[ads_idx], height = z_correc - disp_z,
                                        position = (slab[metal_index].position[0], slab[metal_index].position[1]+disp_y))
                    
                        else:
                            add_adsorbate(slab = slab, adsorbate = adsorbates[ads_idx], height = z_correc - disp_z,
                                        position = (slab[metal_index].position[0], slab[metal_index].position[1]+disp_y2))
                     
                    elif disp_z > 2.5:
                            add_adsorbate(slab = slab, adsorbate = adsorbates[ads_idx], height = z_correc - 2.0,
                                        position = (slab[metal_index].position[0], slab[metal_index].position[1]+disp_y3))
                    
#                    view(slab)
                    for n_atom_ads in range(len(adsorbates[ads_idx])):
                        INCAR['MAGMOM'].append(0.6)
                    
                    mag_dic = dict(zip(slab.get_chemical_symbols(), INCAR['MAGMOM']))
                    slab.get_chemical_symbols
                    INCAR['MAGMOM'] = [mag_dic[x] for x in sorted(slab.get_chemical_symbols())]
                    
                    magm = []
                    for m, g in itertools.groupby(INCAR['MAGMOM'], lambda x: float(x)):
                        magm.append("{}*{}".format(len(tuple(g)),m))
                    
                    bond_length = slab.get_distance(metal_index, 32)
                    
                    f.writelines('%02d_%s_%s\n' % (idx + 1.0, formula, adsorbate_names[ads_idx]))
                    f.writelines('\tformula: %s\n' % (slab.get_chemical_formula()))
                    f.writelines('\tbond length between metal-adsorbate: %1.3f\n' % (bond_length))
                    f.writelines(['\tMAGMOM: ',str(magm),'\n'])
                
                    slab_sorted = sort(slab)
                    del slab[[atom.index > 31 for atom in slab]]
                
                    INCAR['ICHARG'] = 2
                
                    f.writelines('\tLDAUL: %s\tLDAUU: %s\tLDAUJ: %s\n' % (INCAR['LDAUL'],INCAR['LDAUU'],INCAR['LDAUJ']))
                
                    POTCAR = Potcar(potcar_symbols)
                    f.writelines(['\tPOTCAR_symbols: ',str(potcar_symbols), '\n'])
                
#                view(slab_sorted)

                    write_vasp(file_write_path + '%s/POSCAR' % (adsorbate_names[ads_idx]),slab_sorted)
                    INCAR.write_file(file_write_path + '%s/INCAR' % (adsorbate_names[ads_idx]))
                    KPOINTS.write_file(file_write_path + '%s/KPOINTS' % (adsorbate_names[ads_idx]))
                    POTCAR.write_file(file_write_path + '%s/POTCAR' % (adsorbate_names[ads_idx]))
                
                    for n, line in enumerate(fileinput.FileInput('jobscript_vasp.sh')):
                        if '#PBS -N' in line:
                            n_line = n
                        elif '#PBS -q' in line:
                            q_line = n
                        elif '#PBS -l' in line:
                            w_line = n
                
                    PBS_N = '#PBS -N %02d_%s\n' % (idx + 1.0, adsorbate_names[ads_idx])
                    replace_line('jobscript_vasp.sh', n_line, PBS_N)
                
                    queue = 'normal'
                    PBS_q = '#PBS -q %s\n' % (queue)
                    replace_line('jobscript_vasp.sh', q_line, PBS_q)
                
                    walltime = '24:00:00'
                    PBS_w = '#PBS -l walltime=%s\n' % (walltime)
                    replace_line('jobscript_vasp.sh', w_line, PBS_w)
                
                    destination = file_write_path + '/' + adsorbate_names[ads_idx] + '/'
                    job_file = os.getcwd() + '/jobscript_vasp.sh'
                    shutil.copy(job_file, destination)
                
                    f.writelines(['#'*80,'\n'])
                    
                except:
                    f.writelines('%02d_%s_%s result files were not imported\n' % (idx+1.0, formula, adsorbate_names[ads_idx]))
            print('formula %s models are successfully created' % (formula))#v4  
            
        except:
            f.writelines('%02d_%s surface result files were not imported\n' % (idx+1.0, formula))
            f.writelines(['#'*80,'\n'])
            print('failure')
    
        end_time = time.time()
        f.writelines('Excution time for script (sec) : %6.1f\n' % (end_time - start_time))

