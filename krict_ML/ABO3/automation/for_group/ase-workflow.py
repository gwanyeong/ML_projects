# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 15:03:45 2021

@author: gyjung
"""
import time
import os
import shutil

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.io.cif import CifWriter

from ase.calculators.vasp import Vasp
from ase.visualize import view
from ase.io import read
# from ase.dft.kpoints import bandpath    

import matplotlib
matplotlib.use('pdf')

start_time = time.time()

##############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

##############################################################################

MATERIAL_API_KEY = 'your_API_key'
mpr = MPRester(MATERIAL_API_KEY)

TM_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
           'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

alkali_elements = ['Sr', 'Ba', 'K', 'Na', 'Li', 'Rb', 'Cs', 'Be', 'Mg',
                   'Ca', 'Ba', 'Si']


##############################################################################
# Hubbard U choice
U_dict = {'Co':3.32, 'Cr':3.7, 'Fe':5.3,'Mn':3.9, 'Mo':4.38, 'Ni':6.2,
          'V':3.25, 'W':7.17}
U_elements = list(U_dict.keys())

##############################################################################
entries = mpr.query(criteria = {"elements":{"$all":["O"], "$in":TM_list},
                                "anonymous_formula":{"A":1,"B":1,"C":3},
                                "nelements":3,"spacegroup.number":221,
                                "crystal_system":"cubic","spacegroup.symbol":"Pm-3m"},
                    properties = ["materials_id","pretty_formula",
                                  "formation_energy_per_atom","cif","structure",
                                  "band_gap","input.incar","magnetic_type",
                                  "total_magnetization","e_above_hull","volume","theoretical"])

print("%d entries were found" % len(entries))

sorted_entries = sorted(entries, key = lambda e: e['e_above_hull'])
# test_entries = sorted_entries[:10]

##############################################################################

originalPath = os.getcwd()

for idx, entry in enumerate(sorted_entries):
    formula = entry['pretty_formula']

    createFolder(originalPath + '/%02d_%s' % (idx + 1, formula))
    os.chdir(originalPath + '/%02d_%s' % (idx + 1, formula))
    print(os.getcwd())

    struc = entry['structure']
    w = CifWriter(struc)
    w.write_file('%02d_%s.cif' % (idx + 1, formula))

    model = read('%02d_%s.cif' % (idx + 1, formula), primitive_cell = False)
 #   view(model)
    elements = model.get_chemical_symbols()
        
    # MAGMOM settings
    mag_dict ={}
    for el in elements:
        if el in TM_list:
            mag_dict[el] = 5.0  # ferromagnetic
        else:
            mag_dict[el] = 0.6
#   print(mag_dict)
    magmoms = [mag_dict[el] for el in elements] # Edit here if AFM or NM
    
    # Hubbard U setting   
    ldau = False
    ldau_luj = {}
    for el in elements:
        if el in U_elements:
            ldau_luj[el] = {'L':2,'U':U_dict[el],'J':0}
            ldau = True
            
    calc = Vasp(kpts=(6,6,6), system = formula,
                xc = 'pbe', istart = 0, icharg = 1, encut = 520, 
                ediff = 2e-06, lreal = False, algo = 'fast', ediffg = -0.02,
                ismear = -5, sigma = 0.05, lorbit = 11, 
                npar = 16, lplane = True, ncore = 16, ldau = ldau, ldau_luj = ldau_luj)
    
    model.calc = calc
    
    # Spin-restricted SPE
    print("%02d_%s: spin-restricted SPE calc." % (idx + 1, formula))
    createFolder('np')
    calc.set(ispin = 1, nsw = 0, isif = 2, ibrion = -1, directory = 'np')
    model.get_potential_energy()
        
    # Spin-unrestricted Geop
    print("%02d_%s: spin-polar geop." % (idx + 1, formula))
    createFolder('opt')
    shutil.copy('np/CHGCAR', 'opt/') 
    calc.set(ispin = 2, nsw = 200, isif = 3, ibrion = 2, directory = 'opt')
    model.set_initial_magnetic_moments(magmoms = magmoms)
    model.get_potential_energy()
  
    # Spin-unrestricted Geop (2nd)
    print("%02d_%s: spin-polar geop.(2nd)" % (idx + 1, formula))
    createFolder('2nd')
    shutil.copy('opt/CHGCAR', '2nd/')
    calc.set(ispin = 2, nsw = 200, isif = 3, ibrion = 2, directory = '2nd')
    model.get_potential_energy()
    
    # Band structure by SPE.
    print("%02d_%s: band structure calc." % (idx + 1, formula))
    createFolder('bands')
    shutil.copy('2nd/CHGCAR','bands/')
  # kpath = bandpath(path = 'GXMGRXMR',cell = model.cell, npoints = 10)
    kpts = {'path':'GXMGRXMR', 'npoints':60}
    calc.set(nsw = 0, isym = -1, ibrion = -1, isif = 2, directory = 'bands', icharg = 11,
             ismear = 0, nedos = 3000,  kpts = kpts, reciprocal = True)
    model.get_potential_energy()

    # Band structure export
    createFolder('../plots')
    bs = calc.band_structure()
    bs.plot(emin = -10, emax = 10, filename = '../plots/%02d_%s_bs.pdf' % (idx + 1, formula))

end_time = time.time()
print('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
