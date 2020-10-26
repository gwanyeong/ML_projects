# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 01:18:14 2020

@author: gyjung
"""

import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from dotenv import load_dotenv

from pymatgen import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter


###############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)
###############################################################################

load_dotenv('.env')
MATERIAL_API_KEY = os.getenv('MATERIAL_API_KEY')

mpr = MPRester(MATERIAL_API_KEY)

df_entries = pd.read_csv('mpid_list_v3.csv')

###############################################################################

num_ini = 0
num_fin = len(df_entries)

createFolder('results/plots/pourbaix')

with open('results/pbx_analysis.log', 'w') as f:

    for idx in range(num_ini, num_fin):
    
        formula = df_entries['formula'][idx]
        mpid = df_entries['mp_id'][idx]

        print('%03d_%s' % (idx + 1.0, formula))
        f.writelines('%03d_%s\n' % (idx + 1.0, formula))

        try:    
            struc = mpr.get_structure_by_material_id(mpid)
            elements = [str(atom) for atom in struc.species[:2]]
        
            entries = mpr.get_pourbaix_entries(elements)
            entry_ids = [e for e in entries if e.entry_id == mpid]
            if len(entry_ids) == 1:
                pbx = PourbaixDiagram(entries)
                plotter = PourbaixPlotter(pbx)
                entry = entry_ids[0]
                plt = plotter.plot_entry_stability(entry)
                plt.savefig('results/plots/pourbaix/%03d_%s.png' % (idx + 1.0, formula), dpi = 300)
            elif len(entry_ids) > 1:
                f.writelines('%03d_%s: number of entry > 1 \n' % (idx + 1.0, formula))
        except:
            f.writelines('%03d_%s - Error: loading pourbaix entry from MP!\n' % (idx + 1.0, formula))
   
