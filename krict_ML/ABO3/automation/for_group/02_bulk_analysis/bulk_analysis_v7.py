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

# add O2_corr ref. energy
ref_dict['O2_corr'] = -8.45572061

# define U correction values
U_corr_dict = {'V':1.682, 'Cr':2.013, 'Mn':1.68085, 'Fe':2.733,
               'Co':1.874, 'Ni':2.164, 'Mo':3.531, 'W':5.159}

###############################################################################

load_dotenv('.env')
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')

mpr = MPRester(MATERIAL_API_KEY)

df = pd.read_csv('mpid_list_v3.csv')

entries_from_list = mpr.query(criteria = {'material_id':{'$in':list(df.mp_id.values)}},
                              properties = ["material_id","task_id","pretty_formula",
                                            "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                            "structure","input.incar","magnetic_type","total_magnetization",
                                            "e_above_hull","band_gap","volume","theoretical"])

print("Total %d structures were identified from %d mp-ID" % (len(entries_from_list), len(df)))

df_entries_ori = pd.DataFrame(entries_from_list)
df_entries = pd.DataFrame(entries_from_list)


index_list = []
for i in range(len(df)):
    formula = df.formula[i]
    for k, entry in enumerate(df_entries_ori):
        index = df_entries_ori[df_entries_ori['pretty_formula'] == formula].index[0]
        
#   print(i," ",formula," ",index," ",df_entries_ori['pretty_formula'][index])
    index_list.append(index)
    

for n, idx in enumerate(index_list):
    df_entries.iloc[n] = df_entries_ori.iloc[idx]   

"""
atoms_list = []
for i in range(len(df_entries)):
    print(i," ",df_entries['pretty_formula'][i])
    for atom in df_entries['structure'][i].species:
        if str(atom) not in atoms_list:
            atoms_list.append(str(atom))
"""
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
              'VBM', 'E_fermi', 'CBM', 'bgtype', 'bgpath', 'structure'
              ]

Asite_elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']


properties_Asite = ['formula_Na','total_E_Na', 'tot_mag_Na','DHf_Na', 'Ehull_ref_Na', 'Eg_Na', 'Volume_Na',
                    'formula_K','total_E_K', 'tot_mag_K','DHf_K', 'Ehull_ref_K', 'Eg_K', 'Volume_K',
                    'formula_Rb','total_E_Rb', 'tot_mag_Rb','DHf_Rb', 'Ehull_ref_Rb', 'Eg_Rb', 'Volume_Rb',
                    'formula_Cs','total_E_Cs', 'tot_mag_Cs','DHf_Cs', 'Ehull_ref_Cs', 'Eg_Cs', 'Volume_Cs',
                    'formula_Mg','total_E_Mg', 'tot_mag_Mg','DHf_Mg', 'Ehull_ref_Mg', 'Eg_Mg', 'Volume_Mg',
                    'formula_Ca','total_E_Ca', 'tot_mag_Ca','DHf_Ca', 'Ehull_ref_Ca', 'Eg_Ca', 'Volume_Ca',
                    'formula_Sr','total_E_Sr', 'tot_mag_Sr','DHf_Sr', 'Ehull_ref_Sr', 'Eg_Sr', 'Volume_Sr',
                    'formula_Ba','total_E_Ba', 'tot_mag_Ba','DHf_Ba', 'Ehull_ref_Ba', 'Eg_Ba', 'Volume_Ba']

###############################################################################

with open('results/bulk_analysis.log', 'w') as f:
    start_time = time.time()
    f.writelines('No\tformula\tEg_MP(eV)\tEg(eV)\tVBM\tE_fermi\tCBM\tbgtype\tBandgap_path\n')
    
    # insert the number of index(models) for analysis
    num_ini = 0
    num = len(df_entries)    # len(df_entries) 
    
    df_bulk = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                           index = [m for m in range(1,num+1)],
                           columns = properties)
    
    df_Asite = pd.DataFrame(np.array([['NaN' for i in range(len(properties_Asite))] for j in range(num)]),
                            index = [m for m in range(1, num+1)],
                            columns = properties_Asite)

    for idx in range(num_ini, num):
        """
        Import from MP_INCAR
        """
        formula = df_entries['pretty_formula'][idx]
        incar = df_entries['input.incar'][idx]

        df_bulk.MP_ID[idx+1] = df_entries['material_id'][idx]
        df_bulk.formula[idx+1] = formula
        
        df_bulk.PREC[idx+1] = incar['PREC']
        df_bulk.ALGO[idx+1] = incar['ALGO']
        
        try:
            df_bulk.IMIX[idx+1] = incar['IMIX']
        except KeyError:
            df_bulk.IMIX[idx+1] = 4  # default value
            
        df_bulk.ISPIN[idx+1] = incar['ISPIN']
        df_bulk.NELM[idx+1] = incar['NELM']
        df_bulk.IBRION[idx+1] = incar['IBRION']
        
        try:
            df_bulk.EDIFF[idx+1] = incar['EDIFF']
        except KeyError:
            df_bulk.EDIFF[idx+1] = 1e-4 # default value
        
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
        
        """
        Import from MP_results
        """
        struc = df_entries['structure'][idx]
       
        # getting conventional unit cell
        sga = SpacegroupAnalyzer(struc, symprec = 0.1)
        conv_struc = sga.get_conventional_standard_structure()
        if len(conv_struc) != 5:
            print('%s : number of atoms - Error' % file_path)

        scaling = len(conv_struc)/len(struc)
        
        df_bulk.total_E_ref[idx+1] = df_entries['energy'][idx] * scaling
        df_bulk.E_atom[idx+1] = df_entries['energy_per_atom'][idx]
        
        df_bulk.cell_param[idx+1] = struc.lattice.abc
        df_bulk.angle[idx+1] = struc.lattice.angles
        
        try:
            magm_f = struc.site_properties['magmom']
            df_bulk.Magm_ref[idx+1] = magm_f
        except:
            df_bulk.Magm_ref[idx+1] = None
        
        df_bulk.tot_mag_ref[idx+1] = df_entries['total_magnetization'][idx] * scaling
        df_bulk.Mag_O_ref[idx+1] = df_entries['magnetic_type'][idx]
        df_bulk.DHf_ref[idx+1] = df_entries['formation_energy_per_atom'][idx]
        df_bulk.Ehull_ref[idx+1] = df_entries['e_above_hull'][idx]
        df_bulk.Eg_ref[idx+1] = df_entries['band_gap'][idx]
        df_bulk.Theoretical[idx+1] = df_entries['theoretical'][idx]
        
        volume_ref = df_entries['volume'][idx] * scaling
        df_bulk.Volume_ref[idx+1] = volume_ref
                
        # save conventional unit cell
        createFolder('results/bulk_models')        
        conv_struc.to(filename = "results/bulk_models/POSCAR_%s" % (formula))
        conv_struc.to(filename = "results/bulk_models/%s.cif" % (formula))
        
        """
        Import from my results
        """
        if os.path.exists('%03d_%s/2nd/' % (idx + 1.0, formula)):
            if os.path.exists('%03d_%s/2nd/cont/' % (idx + 1.0, formula)):
                file_path = '%03d_%s/2nd/cont/' % (idx + 1.0, formula)
            else:
                file_path = '%03d_%s/2nd/' % (idx + 1.0, formula)
        elif os.path.exists('%03d_%s/2nd_opt/' % (idx + 1.0, formula)):
            if os.path.exists('%03d_%s/2nd_opt/cont/' % (idx + 1.0, formula)):
                file_path = '%03d_%s/2nd_opt/cont/' % (idx + 1.0, formula)
            else:
                file_path = '%03d_%s/2nd_opt/' % (idx + 1.0, formula)
        
        print(file_path)
        
        try:
            oszicar = Oszicar(file_path + 'OSZICAR')
            tot_E = oszicar.ionic_steps[-1]['F']
            tot_mag = oszicar.ionic_steps[-1]['mag']
            
            df_bulk.total_E[idx+1] = tot_E
            df_bulk.tot_mag[idx+1] = tot_mag

            v = Vasprun(file_path + 'vasprun.xml')
            volume = v.as_dict()['output']['crystal']['lattice']['volume']
            
           
            df_bulk.Volume[idx+1] = volume
            df_bulk.Err_V_percent[idx+1] = (volume_ref - volume) / volume_ref * 100

            el = v.as_dict()['elements']
            el.remove('O')

            fE = tot_E - ref_dict[el[0]] - ref_dict[el[1]] - 1.5*ref_dict['O2_corr']
    
            if el[0] or el[1] in U_corr_dict:
                for x in range(len(el)):
                    if el[x] in U_corr_dict:
                        fE = fE - U_corr_dict[el[x]]
            fE = fE/len(conv_struc) # eV/atom
                        
            df_bulk.DHf[idx+1] = fE
            df_bulk.Err_DHf[idx+1] = df_bulk.DHf_ref[idx+1] - df_bulk.DHf[idx+1]
            
            outcar = Outcar(file_path + 'OUTCAR')
            
            magm_list = []
            for atom in range(len(conv_struc)):  
                magm_list.append(outcar.magnetization[atom]['tot'])
            df_bulk.Magm[idx + 1] = magm_list
            
            if tot_mag >= 0.5:    # not a physical definition 
                df_bulk.Mag_O[idx+1] = 'FM'
            elif tot_mag < 0.5:
                df_bulk.Mag_O[idx+1] = 'NM'
            
            bulk_opt = read_vasp(file_path + 'CONTCAR')
            elements = []
            for element in bulk_opt.get_chemical_symbols():
                if element not in elements:
                    elements.append(element)

            struc = Structure.from_file(file_path + 'CONTCAR')
            df_bulk.structure[idx+1] = struc
            
            # Asite dataframe
            for Asite_element in Asite_elements:
                if Asite_element in elements:
                    
                    df_Asite['formula_%s' % (Asite_element)][idx + 1]  = formula
                    df_Asite['total_E_%s' % (Asite_element)][idx + 1] = tot_E
                    df_Asite['tot_mag_%s' % (Asite_element)][idx + 1] = tot_mag
                    df_Asite['DHf_%s' % (Asite_element)][idx + 1] = fE
                    df_Asite['Ehull_ref_%s' % (Asite_element)][idx + 1] = df_entries['e_above_hull'][idx]
                    df_Asite['Volume_%s' % (Asite_element)][idx + 1] = volume
                                    
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

            f.writelines("%03d\t%s\t" % (idx + 1.0, formula))
          
            if (type(vbm) and type(cbm)) == float:
                f.writelines("%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t" % (df_entries['band_gap'][idx], bg, vbm, e_fermi, cbm))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
            elif (vbm and cbm) == None:
                f.writelines("%4.3f\tNone\tNone\t%4.3f\tNone\t" % (df_entries['band_gap'][idx], e_fermi))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))

            """
            Plot DOS & Band
            """
            dosv = Vasprun(file_path + 'DOS/vasprun.xml', parse_dos = True)
            cdos = dosv.complete_dos

            bsdosplot = BSDOSPlotter(bs_projection = "elements", dos_projection = "elements",
                                     vb_energy_range = 10, cb_energy_range = 10, egrid_interval = 2,
                                     font = 'DejaVu Sans')
            bsdosplot.get_plot(bs, cdos).savefig('results/plots/%03d_%s.png' % (idx + 1.0, formula), dpi = 300)
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

    df_bulk.to_csv('results/df_bulk_v2.csv')
    df_bulk.to_pickle('results/df_bulk_v2.pkl')
    
    df_Asite.to_csv('results/df_Asite.csv')
    df_Asite.to_pickle('results/df_Asite.pkl')

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))


