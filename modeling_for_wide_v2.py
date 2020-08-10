# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 08:23:06 2020

@author: gyjung
"""

import os
import shutil

from ase.io.xsd import write_xsd
from ase.io.vasp import read_vasp, write_vasp
from pymatgen.io.vasp.inputs import Incar, Kpoints

# from ase.io.xsd import write_xsd
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
       
cell_params = scalar_sum(zwitter_AA.cell)

tot_num = zwitter_AA.get_global_number_of_atoms()

x_data = []
z_data = []
for i in range(tot_num):
    x_data.append(zwitter_AA[i].x)
    z_data.append(zwitter_AA[i].z)

zwitter_temp = zwitter_AA
# view(zwitter_temp)

# count = 0

INCAR = Incar.from_file('INCAR')
KPOINTS = Kpoints.from_file('KPOINTS')

idx_list = []

scaled_positions = zwitter_temp.get_scaled_positions()

for idx in range(tot_num):
    if scaled_positions[idx][1] > 0.35 and scaled_positions[idx][1] < 0.65:
        idx_list.append(idx)
print('number of atoms in top layer = ',len(idx_list))


a1 = zwitter_AA.cell[0][0] / cell_params[0]
a2 = zwitter_AA.cell[0][1] / cell_params[0]
a3 = zwitter_AA.cell[0][2] / cell_params[0]
print('cell vector A: %4.5f %4.5f %4.5f ' % (a1,a2,a3))

b1 = zwitter_AA.cell[1][0] / cell_params[1]
b2 = zwitter_AA.cell[1][1] / cell_params[1]
b3 = zwitter_AA.cell[1][2] / cell_params[1]
print('cell vector B: %4.5f %4.5f %4.5f ' % (b1,b2,b3))

c1 = zwitter_AA.cell[2][0] / cell_params[2]
c2 = zwitter_AA.cell[2][1] / cell_params[2]
c3 = zwitter_AA.cell[2][2] / cell_params[2]
print('cell vector C: %4.5f %4.5f %4.5f ' % (c1,c2,c3))


for a_int in range(0,11):
    for c_int in range(0,11):
        
        createFolder('%02d_%02d_shift' % (a_int, c_int))
        createFolder('xsd')

        for i in range(tot_num):
            if i in idx_list:
#               count = count + 1
#               print(i," ",zwitter_temp[i]," ",count)
                zwitter_temp[i].x = x_data[i] + a_int*cell_params[0]*a1*0.1 - c_int*cell_params[2]*c1*0.1
                zwitter_temp[i].z = z_data[i] + a_int*cell_params[0]*a3*0.1 - c_int*cell_params[2]*c3*0.1
                
        INCAR.write_file('%02d_%02d_shift/INCAR' % (a_int, c_int))
        KPOINTS.write_file('%02d_%02d_shift/KPOINTS' % (a_int, c_int))    
        write_vasp('%02d_%02d_shift/POSCAR' % (a_int, c_int), zwitter_temp)
        
        write_xsd('xsd/%02d_%02d.xsd' % (a_int, c_int), zwitter_temp)
        
        # jobscript copy
        destination = '%02d_%02d_shift' % (a_int, c_int)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)
       
