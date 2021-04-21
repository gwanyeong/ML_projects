# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:36:16 2021

@author: gyjung
"""

import os
import fileinput
import shutil
import pandas as pd

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
def replace_num(filename, ini_num, fin_num):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if 'num_ini = ' in line:
            ini_line = n
        elif 'num_fin = ' in line:
            fin_line = n
    new_ini_line = 'num_ini = %s\n' % ini_num
    replace_line(filename, ini_line, new_ini_line)

    new_fin_line = 'num_fin = %s\n' % fin_num
    replace_line(filename, fin_line, new_fin_line)
    
##############################################################################
def replace_jobscripts(filename, ini_num, fin_num, formula, queue):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if '#PBS -N' in line:
            n_line = n
        elif '#PBS -q' in line:
            q_line = n
        elif '/scratch/x2045a01/' in line:
            j_line = n
            
    PBS_N = '#PBS -N %03d_%03d_slab\n' % (ini_num + 1, fin_num)
    replace_line(filename, n_line, PBS_N)
    
    PBS_q = '#PBS -q %s\n' % (queue)
    replace_line(filename, q_line, PBS_q)

    command = '/scratch/x2045a01/anaconda3/envs/cgcnn/bin/python \
               jobs/slab_%03d_%03d.py >> logs/stdout_%03d_%03d\n' % (ini_num+1, fin_num, ini_num+1, fin_num)
    replace_line(filename, j_line, command)

##############################################################################

num_ini = 280
num_fin = 300
count = 3

df_entries = pd.read_csv('mpid_list_v3.csv')

filename = 'slab_calc.py'

createFolder('logs')
createFolder('jobs')

for idx in range(num_ini, num_fin, count):
    formula = df_entries['formula'][idx]
    n_ini = idx
    if idx + count < num_fin:
        n_fin = idx + count
    else:
        n_fin = num_fin
    replace_num(filename = filename, ini_num = n_ini, fin_num = n_fin)
    shutil.copy(filename, 'jobs/slab_%03d_%03d.py' % (n_ini + 1, n_fin))
    replace_jobscripts(filename = 'jobscript_ase.sh', 
                       ini_num = n_ini, fin_num = n_fin,
                       formula = formula, queue = 'normal')
    
    os.system('qsub jobscript_ase.sh')
