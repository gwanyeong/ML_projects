# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:00:05 2021

@author: gyjung
"""

import os
import fileinput

##############################################################################
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
   
##############################################################################
def replace_elements(filename, elements):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if 'target_elements = ' in line:
            n_line = n
    new_line = 'target_elements = %s\n' % elements
    replace_line(filename,n_line, new_line)
    
##############################################################################
def replace_jobscripts(filename, label, model_type, queue):
    for n, line in enumerate(fileinput.FileInput(filename)):
        if '#PBS -N' in line:
            n_line = n
        elif '#PBS -q' in line:
            q_line = n
        elif '/scratch/x2045a01/' in line:
            j_line = n
            
    PBS_N = '#PBS -N %s_%s\n' % (model_type, label)
    replace_line(filename, n_line, PBS_N)
    
    PBS_q = '#PBS -q %s\n' % (queue)
    replace_line(filename, q_line, PBS_q)
    
    command = '/scratch/x2045a01/anaconda3/envs/cgcnn/bin/python main.py >> logs/stdout_%s\n' % (label)
    replace_line(filename, j_line, command)

##############################################################################

job_list = [['Sc','Ti'],['V'],['Cr'],['Mn'],['Fe'],['Co'],['Ni'],['Cu','Zn']]

for idx, job in enumerate(job_list):
    TM = ''.join(job_list[idx])
    
    replace_elements('main.py',job)
    replace_jobscripts('jobscript_ase.sh',TM, 'svn3','flat')
    
    os.system('qsub jobscript_ase.sh')

