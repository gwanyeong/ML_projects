# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 08:23:06 2020

@author: gyjung
"""

import os
import shutil

from ase.io.vasp import read_vasp, write_vasp
# from ase.io.xsd import write_xsd
from pymatgen.io.vasp.inputs import Incar, Kpoints

# from ase.visualize import view
from math import sqrt


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

zwitter_AA = read_vasp('POSCAR')
# view(zwitter_AA)


def scalar_sum(array):
    vec_list = []
    for i in range(len(array)):
        vec_list.append(sqrt(array[i][0]**2 + array[i][1]**2 + array[i][2]**2))
    return vec_list   
       
cell_params_ori = scalar_sum(zwitter_AA.cell)

tot_num = zwitter_AA.get_global_number_of_atoms()

#x_data = []
y_data = []
#z_data = []

for i in range(tot_num):
 #   x_data.append(zwitter_AA[i].x)
    y_data.append(zwitter_AA[i].y)
  #  z_data.append(zwitter_AA[i].z)
        
zwitter_temp = zwitter_AA
# view(zwitter_temp)

gap_dist = cell_params_ori[1]/2
dec_ratio = 0.05 # 5% of gap distance

INCAR = Incar.from_file('INCAR')
KPOINTS = Kpoints.from_file('KPOINTS')

idx_list = []

scaled_positions = zwitter_temp.get_scaled_positions()

for idx in range(tot_num):
    if scaled_positions[idx][1] > 0.38 and scaled_positions[idx][1] < 0.62:
        idx_list.append(idx)

# decrease gap distance from 0 to 50% with 5% interval
for y_cnt in range(0,11):
    createFolder('%02d_vertical' % (y_cnt))
    
    for i in range(tot_num):
        if i in idx_list:
            zwitter_temp[i].y = y_data[i] - 0.5*y_cnt*gap_dist*dec_ratio
    
    zwitter_temp.cell[1][1] -= gap_dist*dec_ratio
  # view(zwitter_temp)
    
    INCAR.write_file('%02d_vertical/INCAR' % (y_cnt))
    KPOINTS.write_file('%02d_vertical/KPOINTS' % (y_cnt))    
    write_vasp('%02d_vertical/POSCAR' % (y_cnt), zwitter_temp)
  # write_xsd('%02d_vertical.xsd' % y_cnt, zwitter_temp)
        
    # jobscript copy
    destination = '%02d_vertical' % (y_cnt)
    job_file = os.getcwd() + '/jobscript_vasp.sh'
    shutil.copy(job_file, destination)
       

