# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 08:23:06 2020

@author: gyjung
"""

import os
import shutil

from ase.io.xsd import read_xsd
from ase.io.vasp import read_vasp, write_vasp
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
       
cell_params = scalar_sum(zwitter_AA.cell)

tot_num = zwitter_AA.get_global_number_of_atoms()

x_data = []
z_data = []
for i in range(tot_num):
    x_data.append(zwitter_AA[i].x)
    z_data.append(zwitter_AA[i].z)
        
zwitter_temp = zwitter_AA
# view(zwitter_temp)


INCAR = Incar.from_file('INCAR')
KPOINTS = Kpoints.from_file('KPOINTS')

#count = 0
for x_int in range(0,6):
    for z_int in range(0,6):
        createFolder('%d_%d_shift' % (x_int, z_int))

        for i in range(tot_num):
            
            if zwitter_temp[i].y > 4.5 and zwitter_temp[i].y < 7.0:
                #count = count + 1
                #print(i," ",zwitter_temp[i]," ",count)
                zwitter_temp[i].x = x_data[i] + x_int*cell_params[0]*0.01
                zwitter_temp[i].z = z_data[i] + z_int*cell_params[2]*0.01

        INCAR.write_file('%d_%d_shift/INCAR' % (x_int, z_int))
        KPOINTS.write_file('%d_%d_shift/KPOINTS' % (x_int, z_int))    
        write_vasp('%d_%d_shift/POSCAR' % (x_int, z_int), zwitter_temp)
        
        # jobscript copy
        destination = '%d_%d_shift' % (x_int, z_int)
        job_file = os.getcwd() + '/jobscript_vasp.sh'
        shutil.copy(job_file, destination)
       



