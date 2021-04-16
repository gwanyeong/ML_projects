# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:00:05 2021

@author: gyjung
"""

import os
import fileinput
import shutil

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
def replace_model(filename, elements, model_type):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if 'target_elements = ' in line:
            n_line = n
        elif 'model_type = ' in line:
            m_line = n
    new_el_line = 'target_elements = %s\n' % elements
    replace_line(filename, n_line, new_el_line)

    new_model_line = "model_type = '%s'\n" % model_type
    replace_line(filename, m_line, new_model_line)
    
##############################################################################
def replace_jobscripts(filename, label, model_type, queue):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if '#PBS -N' in line:
            n_line = n
        elif '#PBS -q' in line:
            q_line = n
        elif '/scratch/x2045a01/' in line:
            j_line = n
            
    PBS_N = '#PBS -N %s_%s_cnt\n' % (model_type, label)
    replace_line(filename, n_line, PBS_N)
    
    PBS_q = '#PBS -q %s\n' % (queue)
    replace_line(filename, q_line, PBS_q)

    command = '/scratch/x2045a01/anaconda3/envs/cgcnn/bin/python jobs/calc_%s_cnt.py >> logs/stdout_%s_cnt\n' % (label, label)
    replace_line(filename, j_line, command)

##############################################################################

# 3d TMs
# job_list = [['Sc','Ti'],['V'],['Cr'],['Mn'],['Fe'],['Co'],['Ni'],['Cu','Zn']]

# 4d TMs
job_list = [['Y','Zr'],['Nb'],['Mo'],['Tc'],['Ru'],['Rh'],['Pd'],['Ag','Cd']]

# 5d TMs
# job_list = [['Lu','Hf'],['Ta'],['W'],['Re'],['Os'],['Ir'],['Pt'],['Au','Hg']]


# filename = 'main.py'
filename = 'calculate_cont.py'

createFolder('logs')
createFolder('jobs')
for idx, els in enumerate(job_list):
    TM = ''.join(job_list[idx])
    
    replace_model(filename = filename, elements = els, model_type = 'dvn4')
    shutil.copy(filename, 'jobs/calc_%s_cnt.py' % TM)
    replace_jobscripts(filename = 'jobscript_ase.sh', label = TM, model_type = 'dvn4', queue = 'flat')
    
    os.system('qsub jobscript_ase.sh')
