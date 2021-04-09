# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 21:18:38 2021

@author: gyjung
"""

import os

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

originalPath = os.getcwd()
Gn_s[24].symbol = 'Fe'
view(Gn_s)
write(originalPath + '/models/svc3_ini.cif', images = Gn_s)

b_ligands = [11, 22, 26]

# SVN1C2
for atom in Gn_s:
    if atom.index == 11:
        atom.symbol = 'B'
view(Gn_s)
write(originalPath + '/models/svb1c2_ini.cif', images = Gn_s)

# SVN2C1        
for atom in Gn_s:
    if atom.index == 22:
        atom.symbol = 'B'
view(Gn_s)
write(originalPath + '/models/svb2c1_ini.cif', images = Gn_s)

# SVN3        
for atom in Gn_s:
    if atom.index == 26:
        atom.symbol = 'B'
view(Gn_s)
write(originalPath + '/models/svb3_ini.cif', images = Gn_s)

target_elements = ['Sc']

# Modeling
for idx, TM in enumerate(TM_elements):
    if TM in target_elements:
    
        formula = TM + '_svb3'
        model = read(originalPath + '/models/svb3_ini.cif')
        model[24].symbol = TM

        model = sort(model)
        write(originalPath + '/models/%02d_%s_svb3.cif' % (idx + 1, TM), images = model)

        view(model)


