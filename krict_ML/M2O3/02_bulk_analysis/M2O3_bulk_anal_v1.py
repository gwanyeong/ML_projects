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

from dotenv import load_dotenv

from ase.io.vasp import read_vasp

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

# define U correction values
U_corr_dict = {'V':1.682, 'Cr':2.013, 'Mn':1.68085, 'Fe':2.733,
               'Co':1.874, 'Ni':2.164, 'Mo':3.531, 'W':5.159}

###############################################################################

load_dotenv('.env')
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')
mpr = MPRester(MATERIAL_API_KEY)

TM_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
           'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']

df = pd.read_csv('mpid_list_v2.csv')

M2O3_list = []
for TM in TM_list:
    M2O3 = TM + "2O3"
    M2O3_list.append(M2O3)
    
entries = mpr.query(criteria = {'material_id':{'$in':list(df.mp_id.values)}},
                              properties = ["material_id","task_id","pretty_formula",
                                            "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                            "structure","input.incar","magnetic_type","total_magnetization",
                                            "e_above_hull","band_gap","volume","theoretical"])

print("Total %d structures were identified from %d mp-ID" % (len(entries), len(df)))

entries = sorted(entries, key = lambda e: ['e_above_hull'])

df_entries = pd.DataFrame(entries)

df_entries = df_entries.drop_duplicates(['pretty_formula'], keep = 'first')
print('%d species are found among %d entries' % (len(df_entries),len(M2O3_list)))

MP_list = df_entries['pretty_formula'].tolist() # [['Sc2O3', 'Y2O3', 'Fe2O3', 'Rh2O3', 'V2O3', 'Ti3O2', 'Mn2O3', 'Ti2O3', 'Cr2O3']]
print(MP_list)


MP_dict = {key:value for value, key in enumerate(MP_list)}
print(MP_dict)
###############################################################################

createFolder('results')
createFolder('results/plots')

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
    #print(df_bulk)
    
    for idx, formula in enumerate(M2O3_list):
        if formula in MP_list:            
            #import MP_INCAR
            
            idx2 = MP_dict[formula]
            formula2 = df_entries['pretty_formula'][idx2]
                         
            print('index in df_bulk(%s): '%(formula), idx+1, 'index in MP_list(%s): '%(formula2),idx2)
            incar = df_entries['input.incar'][idx2]
                    
            df_bulk.MP_ID[idx+1] = df_entries['material_id'][idx2]
            df_bulk.formula[idx+1] = formula2
        
            df_bulk.PREC[idx+1] = incar['PREC']
            df_bulk.ALGO[idx+1] = incar['ALGO']
            
            try:
                df_bulk.IMIX[idx+1] = incar['IMIX']
            except KeyError:
                df_bulk.IMIX[idx+1] = 4
            
            try:
                df_bulk.ISPIN[idx+1] = incar['ISPIN']
            except:
                df_bulk.ISPIN[idx+1] = 3
                
            df_bulk.NELM[idx+1] = incar['NELM']
            df_bulk.IBRION[idx+1] = incar['IBRION']
            
            try:
                df_bulk.EDIFF[idx+1] = incar['EDIFF']
            except:
                df_bulk.EDIFF[idx+1] = 1e-4
            
            df_bulk.ISIF[idx+1] = incar['ISIF']
            df_bulk.NSW[idx+1] = incar['NSW']
            df_bulk.ENCUT[idx+1] = incar['ENCUT']
            df_bulk.MAGMOM[idx+1] = incar['MAGMOM']
            df_bulk.ISMEAR[idx+1] = incar['ISMEAR']
            df_bulk.SIGMA[idx+1] = incar['SIGMA']
            
            
            
            try:
                df_bulk.LDAU[idx+1] = incar['LDAU']
                df_bulk.LDAUU[idx+1] = incar['LDAUU']
            except KeyError:
                df_bulk.LDAU[idx+1] = None
                df_bulk.LDAUU[idx+1] = None    
            
            
            #import MP_results
            
            struc_MP = df_entries['structure'][idx2]       
            sga = SpacegroupAnalyzer(struc_MP, symprec = 0.1)
            conv_struc = sga.get_conventional_standard_structure()
            if len(conv_struc) != 30:
                 print('%s : number of atoms - Error' % idx)

            scaling = len(conv_struc)/len(struc_MP)
            print('scaling factor for %s: ' % (formula),scaling)
        
            df_bulk.total_E_ref[idx+1] = df_entries['energy'][idx2] * scaling
            df_bulk.E_atom[idx+1] = df_entries['energy_per_atom'][idx2]
        
            df_bulk.cell_param[idx+1] = struc_MP.lattice.abc
            df_bulk.angle[idx+1] = struc_MP.lattice.angles
        
            try:
                magm_f = struc_MP.site_properties['magmom']
                df_bulk.Magm_ref[idx+1] = magm_f
            except:
                df_bulk.Magm_ref[idx+1] = None
        
            df_bulk.tot_mag_ref[idx+1] = df_entries['total_magnetization'][idx2] * scaling
            df_bulk.Mag_O_ref[idx+1] = df_entries['magnetic_type'][idx2]
            df_bulk.DHf_ref[idx+1] = df_entries['formation_energy_per_atom'][idx2]
            df_bulk.Ehull_ref[idx+1] = df_entries['e_above_hull'][idx2]
            df_bulk.Eg_ref[idx+1] = df_entries['band_gap'][idx2]
            df_bulk.Theoretical[idx+1] = df_entries['theoretical'][idx2]
        
            volume_ref = df_entries['volume'][idx2] * scaling
            df_bulk.Volume_ref[idx+1] = volume_ref
                
            createFolder('results/bulk_models')        
            conv_struc.to(filename = "results/bulk_models/POSCAR_%s" % (formula2))
            conv_struc.to(filename = "results/bulk_models/%s.cif" % (formula2))
            
            print(df_bulk)
            
            
        #import my results
        file_path = '%02d_%s/cont/' % (idx + 1.0, formula)
        print(file_path)
        
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
            if formula in MP_list:
                df_bulk.Err_V_percent[idx+1] = (volume_ref - volume) / volume_ref * 100
            bulk_opt = read_vasp(file_path + 'CONTCAR')
            
            elements = []
            for element in bulk_opt.get_chemical_symbols():
                if element not in elements:
                    elements.append(element)

            struc = Structure.from_file(file_path + 'CONTCAR')
            df_bulk.structure[idx+1] = bulk_opt
            
            el = v.as_dict()['elements']
            el.remove('O')

            fE = tot_E - 6*(2*ref_dict[el[0]] + 1.5*ref_dict['O2_corr'])
    
            if el[0] in U_corr_dict:
               fE = fE - 6*(2*U_corr_dict[el[0]])
            fE = fE/len(struc) # eV/atom
                        
            df_bulk.DHf[idx+1] = fE
            if formula in MP_list:
                df_bulk.Err_DHf[idx+1] = df_bulk.DHf_ref[idx+1] - df_bulk.DHf[idx+1]
                        
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
            # print('%03d_%s/2nd files are not found' % (idx + 1.0, formula))
            f.writelines('%s files are not found\n' % (file_path))
            continue

        try:
            bsv = BSVasprun(file_path + 'DOS/vasprun.xml', parse_projected_eigen = True)
            bs = bsv.get_band_structure(kpoints_filename = file_path + 'DOS/KPOINTS', line_mode = True)
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
            dosv = Vasprun(file_path + 'DOS/vasprun.xml', parse_dos = True)
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




