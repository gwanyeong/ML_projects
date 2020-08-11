# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 08:32:56 2020

@author: gyjung
"""

# Analysis for VASP results
import csv

from pymatgen import MPRester

from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.io.vasp.outputs import Oszicar

import pandas as pd

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

mpr = MPRester('mwEnJkSB9Cusgdf1')
mpid_list = []

with open('mpid_list.csv', 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for line in reader:
        mpid = line[0]
        mpid_list.append(mpid)

entries = mpr.query(criteria = {'material_id':{'$in':mpid_list}},
                    properties = ['material_id', 'pretty_formula', 'band_gap', 'e_above_hull'])

print(len(entries)," structures were identified among ", len(mpid_list))

sorted_entries = sorted(entries, key = lambda e: e['e_above_hull'])

with open('summary', 'w') as f:
    f.writelines('No\tformula\ttotal E\tMagm\tF.E.(eV/f.u.)\tVolume\t')
    f.writelines('VBM\tE_fermi\tCBM\tbgtype\tBand_gap(eV)\tBand_gap_MP(eV)\tBandgap_path\n')

    for idx in range(len(sorted_entries)):
        formula = sorted_entries[idx]['pretty_formula']
        bg_ref = sorted_entries[idx]['band_gap']
    
        try:
            oszicar = Oszicar('%03d_%s/2nd/OSZICAR' % (idx + 1.0, formula))
            tot_E = oszicar.ionic_steps[-1]['F']
            tot_mag = oszicar.ionic_steps[-1]['mag']

            v = Vasprun('%03d_%s/2nd/vasprun.xml' % (idx + 1.0, formula))
            volume = v.as_dict()['output']['crystal']['lattice']['volume']

            el = v.as_dict()['elements']
            el.remove('O')

            fE = tot_E - ref_dict[el[0]] - ref_dict[el[1]] - 1.5*ref_dict['O2_corr']
    
            if el[0] or el[1] in U_corr_dict:
                for x in range(len(el)):
                    if el[x] in U_corr_dict:
                        fE = fE - U_corr_dict[el[x]]

            #print(i+1," ", formula," ",tot_E," ",tot_mag," ",fE," ",volume) 
            f.writelines('%d\t%s\t%4.4f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t' % (idx+1.0, formula, tot_E, tot_mag, fE, volume, bg_ref))

        except FileNotFoundError:
            print('%03d_%s/2nd files are not found' % (idx + 1.0, formula))
            f.writelines('%03d_%s/2nd files are not found\n' % (idx + 1.0, formula))
            continue

        try:
            bsv = BSVasprun('%03d_%s/2nd/DOS/vasprun.xml' % (idx+1.0, formula), parse_projected_eigen = True)
            bs = bsv.get_band_structure(kpoints_filename = '%03d_%s/2nd/DOS/KPOINTS' % (idx + 1.0, formula), line_mode = True)
            vbm = bs.as_dict()['vbm']['energy']
            efermi = bs.as_dict()['efermi']
            cbm = bs.as_dict()['cbm']['energy']
 
            if bs.get_band_gap()['direct'] == True:
                bgtype = 'Direct'
            elif bs.get_band_gap()['direct'] == False:
                bgtype = 'Indirect'
 
            bg = bs.get_band_gap()['energy']
            bgpath = bs.get_band_gap()['transition']
          
            if (type(vbm) and type(cbm)) == float:
                f.writelines("%4.3f\t%4.3f\t%4.3f\t%4.3f\t" % (bg, vbm, efermi, cbm))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
            elif (vbm and cbm) == None:
                f.writelines("None\tNone\t%4.3f\tNone\t" % (efermi))
                f.writelines("%s\t%s\n" % (bgtype, bgpath))
    
        except FileNotFoundError:
            print('%03d_%s/2nd/DOS files are not found' % (idx + 1.0, formula))
            f.writelines('-\t-\t-\t-\t-\t-\n')
            continue

        except:
            print('%03d_%s/2nd/DOS files - unknown error' % (idx + 1.0, formula))
            f.writelines('err\terr\terr\terr\terr\terr\n')
            continue

