#!/usr/bin/env python
# coding: utf-8

# In[23]:


# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 00:44:34 2020
@author: gyjung
"""

import os
import csv
import pickle
import time
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

from ase.io.vasp import read_vasp
from ase.io.xsd import write_xsd

from pymatgen import MPRester, Structure
from pymatgen.io.vasp import Vasprun, BSVasprun, Oszicar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.plotter import BSDOSPlotter

###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

###############################################################################
with open('form_E_ref.csv', 'r') as f:
    ref_data = pd.read_csv(f)
    ref_dict = dict.fromkeys(ref_data['element'])
    for i in range(len(ref_data)):
        ref_dict.update({ref_data['element'][i]:ref_data['E'][i]})

# add O2_corr, H2_corr ref. energy
ref_dict['O2_corr'] = -8.45572061
ref_dict['H2_corr'] = -6.781

# define U correction values
U_corr_dict = {'V':1.682, 'Cr':2.013, 'Mn':1.68085, 'Fe':2.733,
               'Co':1.874, 'Ni':2.164, 'Mo':3.531, 'W':5.159}

###############################################################################

TM_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
           'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

MHO2_list = []
for TM in TM_list:
    MHO2 = TM + "HO2"
    MHO2_list.append(MHO2)

###############################################################################
createFolder('results')
createFolder('results/plots')
createFolder('results/bulk_models')        

properties = ['MP_ID', 'formula',
              'PREC', 'ALGO', 'ISPIN', 'IMIX', 'NELM', 'IBRION', 'EDIFF', 'NSW',
              'ISIF', 'ENCUT', 'MAGMOM', 'ISMEAR', 'SIGMA', 'LDAU', 'LDAUU',
              'total_E_ref', 'E_atom', 'total_E', 'cell_param', 'angle',
              'Magm_ref', 'Mag_O_ref','Magm','Mag_O','tot_mag_ref','tot_mag',
              'DHf_ref', 'DHf', 'Err_DHf', 'Ehull_ref', 'Eg_ref', 'Eg',
              'Volume_ref', 'Volume', 'Err_V_percent', 'Theoretical',
              'VBM', 'E_fermi', 'CBM', 'bgtype', 'bgpath', 'structure']

###############################################################################
with open('results/bulk_analysis.log','w') as f:
    start_time = time.time()
    f.writelines('No\tformula\tEg_(eV)\tVBM\tE_fermi\tCBM\tbgtype\tBandgap_path\n')
    
    num_ini = 0
    num = len(TM_list)
    
    df_bulk = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                          index = [m for m in range(1, num+1)],columns = properties)
#    print(df_bulk)
    
    for idx, formula in enumerate(MHO2_list):
        
        #import my results
        file_path = '%02d_%s/cont/' % (idx + 1.0, formula)
        
        try:
            oszicar = Oszicar(file_path + 'OSZICAR')
            tot_E = oszicar.ionic_steps[-1]['F']
            tot_mag = oszicar.ionic_steps[-1]['mag']
            
            df_bulk.formula[idx+1] = formula
            df_bulk.total_E[idx+1] = tot_E
            df_bulk.tot_mag[idx+1] = tot_mag

            v = Vasprun(file_path + 'vasprun.xml')
            volume = v.as_dict()['output']['crystal']['lattice']['volume']

            df_bulk.Volume[idx+1] = volume
#            if formula in MP_list:
#                df_bulk.Err_V_percent[idx+1] = (volume_ref - volume) / volume_ref * 100
            bulk_opt = read_vasp(file_path + 'CONTCAR')

            write_xsd('results/bulk_models/poscar_%02d_%s.xsd' % (idx +  1.0, formula), bulk_opt)
    	
            bulk_opt = read_vasp(file_path + 'CONTCAR')
            elements = []
            df_bulk.Volume[idx+1] = volume
            for element in bulk_opt.get_chemical_symbols():
                if element not in elements:
                    elements.append(element)

            struc = Structure.from_file(file_path + 'CONTCAR')
            df_bulk.structure[idx+1] = bulk_opt
            
            el = v.as_dict()['elements']
            el.remove('O')
            el.remove('H')

            fE = tot_E - 2*(ref_dict[el[0]] + 0.5*ref_dict['H2_corr'] + ref_dict['O2_corr'])
    
            if el[0] in U_corr_dict:
                fE = fE - 2*(U_corr_dict[el[0]])
            fE = fE/len(struc) # eV/atom
                        
            df_bulk.DHf[idx+1] = fE
#            if formula in MP_list:
#               df_bulk.Err_DHf[idx+1] = df_bulk.DHf_ref[idx+1] - df_bulk.DHf[idx+1]
                        
            outcar = Outcar(file_path + 'OUTCAR')
            
            magm_list = []
            for atom in range(len(struc)):  
                magm_list.append(outcar.magnetization[atom]['tot'])
            df_bulk.Magm[idx + 1] = magm_list
            
            if tot_mag >= 0.5:    # not a physical definition 
                df_bulk.Mag_O[idx+1] = 'FM'
            elif tot_mag < 0.5:
                df_bulk.Mag_O[idx+1] = 'NM'
                                                
        except FileNotFoundError:
            # print('%02d_%s/cont files are not found' % (idx + 1.0, formula))
            f.writelines('%s files are not found\n' % (file_path))
            continue

        try:
            bsv = BSVasprun(file_path + 'SPE/DOS/vasprun.xml', parse_projected_eigen = True)
            bs = bsv.get_band_structure(kpoints_filename = file_path + 'SPE/DOS/KPOINTS', line_mode = True)
            vbm = bs.as_dict()['vbm']['energy']
            e_fermi = bs.as_dict()['efermi']
            cbm = bs.as_dict()['cbm']['energy']
            
            df_bulk.VBM[idx+1] = vbm
            df_bulk.E_fermi[idx+1] = e_fermi
            df_bulk.CBM[idx+1] = cbm

            bg = bs.get_band_gap()['energy']
            bgpath = bs.get_band_gap()['transition']
            
            if bs.get_band_gap()['direct'] == True:
                bgtype = 'Direct'
            elif bs.get_band_gap()['direct'] == False:
                bgtype = 'Indirect'
            elif bg is None:
                bg = 0.0
                bgtype = None
            
            df_bulk.Eg[idx+1] = bg
            df_bulk.bgpath[idx+1] = bgpath
            df_bulk.bgtype[idx+1] = bgtype

            f.writelines("%02d\t%s\t" % (idx + 1.0, formula))

            if (type(vbm) and type(cbm)) == float:
                f.writelines("%4.3f\t%4.3f\t%4.3f\t%4.3f\t" % (bg, vbm, e_fermi, cbm))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
            elif (vbm and cbm) == None:
                f.writelines("None\tNone\t%4.3f\tNone\t" % (e_fermi))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))

            #Plot DOS & Band
            dosv = Vasprun(file_path + 'SPE/DOS/vasprun.xml', parse_dos = True)
            cdos = dosv.complete_dos

            bsdosplot = BSDOSPlotter(bs_projection = "elements", dos_projection = "elements",
                                     vb_energy_range = 10, cb_energy_range = 10, egrid_interval = 2,
                                     font = 'DejaVu Sans')
            bsdosplot.get_plot(bs, cdos).savefig('results/plots/%02d_%s.png' % (idx + 1.0, formula), dpi = 300)
            if idx < 20:
                print(df_bulk.head(20))
            else:
                print(df_bulk[idx-20:idx])

        except FileNotFoundError:
            f.writelines('%s are not found\n' % (file_path))
            print('%s are not found' % (file_path))
            continue
        
        except:
            f.writelines('%sDOS files - unknown error\n' % (file_path))
            print('%sDOS files - unknown error' % (file_path))
            continue

    df_bulk.to_csv('results/df_bulk.csv')
    df_bulk.to_pickle('results/df_bulk.pkl')

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
