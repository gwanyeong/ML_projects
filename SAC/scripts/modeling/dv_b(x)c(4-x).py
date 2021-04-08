# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 19:59:03 2021

@author: gyjung
"""

import os
# import numpy as np

from ase import Atom
from ase.io import read, write
from ase.visualize import view
from ase.build import make_supercell, sort


##############################################################################
TM_elements = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
               'La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

##############################################################################

Gn = read('models/Gn.cif')
Gn_s = make_supercell(Gn, [[2,0,0],[0,4,0],[0,0,1]])

positions = Gn_s.get_positions()
for pos in positions:
    pos[2] = 0

Gn_s.set_positions(positions)
# view(Gn_super)

originalPath = os.getcwd()

position = (Gn_s[11].position + Gn_s[24].position)/2

del_list = [11,24]
del Gn_s[[atom.index for atom in Gn_s if atom.index in del_list]]
Gn_s.append(Atom('Fe', position= position))

# DVC4
# view(Gn_s)
write(originalPath + '/models/dvc4_ini.cif', images = Gn_s)

n_ligands = [5,9,21,24]

# DVC3N1
for atom in Gn_s:
    if atom.index == 5:
        atom.symbol = 'B'
# view(Gn_s)
write(originalPath + '/models/dvb1c3_ini.cif', images = Gn_s)

# DVC2N2        
for atom in Gn_s:
    if atom.index == 24:
        atom.symbol = 'B'
# view(Gn_s)
write(originalPath + '/models/dvb2c2_ini.cif', images = Gn_s)

# DVC1N3                
for atom in Gn_s:
    if atom.index == 21:
        atom.symbol = 'B'
# view(Gn_s)
write(originalPath + '/models/dvb3c1_ini.cif', images = Gn_s)

# DVN4
for atom in Gn_s:
    if atom.index == 9:
        atom.symbol = 'B'
# view(Gn_s)
write(originalPath + '/models/dvb4_ini.cif', images = Gn_s)

target_elements = ['Sc']

# Modeling
for idx, TM in enumerate(TM_elements):
    if TM in target_elements:
    
        formula = TM + '_dvb4'
        model = read(originalPath + '/models/dvb4_ini.cif')
        model[-1].symbol = TM

        model = sort(model)
        write(originalPath + '/models/%02d_%s_dvb4.cif' % (idx + 1, TM), images = model)

  #     view(model)

