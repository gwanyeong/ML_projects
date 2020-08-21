"""
Created on Tue Aug 20 10:25:52 2020

@author: gyjung
"""

import os
import csv
import pickle
import time

import pandas as pd
import numpy as np

from dotenv import load_dotenv

from pymatgen import MPRester
from pymatgen.io.vasp import Vasprun, BSVasprun, Oszicar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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

mpid_list = []
with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

entries_from_list = mpr.query(criteria = {'material_id':{'$in':mpid_list}},
                              properties = ["material_id","task_id","pretty_formula",
                                            "formation_energy_per_atom","cif", "energy","energy_per_atom",
                                            "structure","band_gap","input.incar","magnetic_type","total_magnetization",
                                            "e_above_hull","band_gap","volume","theoretical"])

print("Total %d structures were identified from %d mp-ID" % (len(entries_from_list), len(mpid_list)))

entries = sorted(entries_from_list, key = lambda e: e['e_above_hull'])

###############################################################################

createFolder('results')

properties = ['MP_ID', 'formula',
              'PREC', 'ALGO', 'ISPIN', 'IMIX', 'NELM', 'IBRION', 'EDIFF', 'NSW',
              'ISIF', 'ENCUT', 'MAGMOM', 'ISMEAR', 'SIGMA', 'LDAU', 'LDAUU',
              'total_E_ref', 'E_atom', 'total_E', 'cell_param', 'angle',
              'Magm_ref', 'Mag_O_ref','Magm','Mag_O','tot_mag_ref','tot_mag',
              'DHf_ref', 'DHf', 'Err_DHf', 'Ehull_ref', 'Eg_ref', 'Eg',
              'Volume_ref', 'Volume', 'Err_V_percent', 'Theoretical',
              'VBM', 'E_fermi', 'CBM', 'bgtype', 'bgpath'
              ]

###############################################################################


with open('results/bulk_analysis.log', 'w') as f:
    start_time = time.time()
    f.writelines('No\tformula\ttotal E\tMagm\tF.E.(eV/f.u.)\tVolume\t')
    f.writelines('Eg_MP(eV)\tEg(eV)\tVBM\tE_fermi\tCBM\tbgtype\tBandgap_path\n')
    
    # insert the number of index(models) for analysis
    num = len(entries)    # len(entries) 
    
    df_bulk = pd.DataFrame(np.array([['NaN' for i in range(len(properties))] for j in range(num)]),
                           index = [m for m in range(1,num+1)],
                           columns = properties)

    for idx in range(num):
        """
        Import from MP_INCAR
        """
        formula = entries[idx]['pretty_formula']
        incar = entries[idx]['input.incar']
        
        df_bulk.MP_ID[idx+1] = entries[idx]['material_id']
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
        struc = entries[idx]['structure']
       
        # getting conventional unit cell
        sga = SpacegroupAnalyzer(struc)
        conv_struc = sga.get_conventional_standard_structure()

        scaling = len(conv_struc)/len(struc)
        
        df_bulk.total_E_ref[idx+1] = entries[idx]['energy'] * scaling
        df_bulk.E_atom[idx+1] = entries[idx]['energy_per_atom']
        
        df_bulk.cell_param[idx+1] = struc.lattice.abc
        df_bulk.angle[idx+1] = struc.lattice.angles
        
        try:
            magm_f = struc.site_properties['magmom']
            df_bulk.Magm_ref[idx+1] = magm_f
        except:
            df_bulk.Magm_ref[idx+1] = None
        
        df_bulk.tot_mag_ref[idx+1] = entries[idx]['total_magnetization'] * scaling
        df_bulk.Mag_O_ref[idx+1] = entries[idx]['magnetic_type']
        df_bulk.DHf_ref[idx+1] = entries[idx]['formation_energy_per_atom']
        df_bulk.Ehull_ref[idx+1] = entries[idx]['e_above_hull']
        df_bulk.Eg_ref[idx+1] = entries[idx]['band_gap']
        df_bulk.Theoretical[idx+1] = entries[idx]['theoretical']
        
        volume_ref = entries[idx]['volume'] * scaling
        df_bulk.Volume_ref[idx+1] = volume_ref
                
        # save conventional unit cell
        createFolder('results/bulk_models')        
        conv_struc.to(filename = "results/bulk_models/POSCAR_%s" % (formula))
        conv_struc.to(filename = "results/bulk_models/%s.cif" % (formula))
        
        """
        Import from my results
        """
        try:
            oszicar = Oszicar('%03d_%s/2nd/OSZICAR' % (idx + 1.0, formula))
            tot_E = oszicar.ionic_steps[-1]['F']
            tot_mag = oszicar.ionic_steps[-1]['mag']
            
            df_bulk.total_E[idx+1] = tot_E
            df_bulk.tot_mag[idx+1] = tot_mag

            v = Vasprun('%03d_%s/2nd/vasprun.xml' % (idx + 1.0, formula))
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
                        
            df_bulk.DHf[idx+1] = fE/5.0
            df_bulk.Err_DHf[idx+1] = df_bulk.DHf_ref[idx+1] - df_bulk.DHf[idx+1]
            
            outcar = Outcar('%03d_%s/2nd/OUTCAR' % (idx + 1.0, formula))
            
            magm_list = []
            for atom in range(len(conv_struc)):  
                magm_list.append(outcar.magnetization[atom]['tot'])
            df_bulk.Magm[idx + 1] = magm_list
            
            if tot_mag >= 0.5:    # not a physical definition 
                df_bulk.Mag_O[idx+1] = 'FM'
            elif tot_mag < 0.5:
                df_bulk.Mag_O[idx+1] = 'NM'
                
        except FileNotFoundError:
            # print('%03d_%s/2nd files are not found' % (idx + 1.0, formula))
            f.writelines('%03d_%s/2nd files are not found\n' % (idx + 1.0, formula))
            continue

        try:
            bsv = BSVasprun('%03d_%s/2nd/DOS/vasprun.xml' % (idx+1.0, formula), parse_projected_eigen = True)
            bs = bsv.get_band_structure(kpoints_filename = '%03d_%s/2nd/DOS/KPOINTS' % (idx + 1.0, formula), line_mode = True)
            vbm = bs.as_dict()['vbm']['energy']
            e_fermi = bs.as_dict()['efermi']
            cbm = bs.as_dict()['cbm']['energy']
            
            df_bulk.VBM[idx+1] = vbm
            df_bulk.E_fermi[idx+1] = e_fermi
            df_bulk.CBM[idx+1] = cbm
            
 
            if bs.get_band_gap()['direct'] == True:
                bgtype = 'Direct'
            elif bs.get_band_gap()['direct'] == False:
                bgtype = 'Indirect'
            
            df_bulk.bgtype[idx+1] = bgtype
 
            bg = bs.get_band_gap()['energy']
            bgpath = bs.get_band_gap()['transition']
            
            df_bulk.bgtype[idx+1] = bgtype
            df_bulk.Eg[idx+1] = bg
            df_bulk.bgpath[idx+1] = bgpath
          
            if (type(vbm) and type(cbm)) == float:
                f.writelines("%4.3f\t%4.3f\t%4.3f\t%4.3f\t" % (bg, vbm, e_fermi, cbm))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
            elif (vbm and cbm) == None:
                f.writelines("None\tNone\t%4.3f\tNone\t" % (e_fermi))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
    
        except FileNotFoundError:
            f.writelines('%03d_%s/2nd/DOS files are not found\n' % (idx + 1.0, formula))
            continue
        
        except:
            f.writelines('%03d_%s/2nd/DOS files - unknown error\n' % (idx + 1.0, formula))
            continue
        
    df_bulk.to_csv('results/df_bulk.csv')
    df_bulk.to_pickle('results/df_bulk.pkl')

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
