# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 14:53:07 2020

@author: gyjung
"""

import os, time, requests
import pandas as pd
from dotenv import load_dotenv

load_dotenv(".env")

### URL and Credential
KRICTDB_URL = "http://calc.chemdx.org"
API_KEY = os.getenv('KRICT_API_KEY')
API_PASS = os.getenv('KRICT_API_PASS')
####
access = "private"    # personal authorization
####

auth = {"KRICTDB_URL": KRICTDB_URL, "API_KEY": API_KEY, "API_PASS": API_PASS}
data = {"api_key": auth['API_KEY'], "api_pass": auth['API_PASS'],
        "access": access, "groupid": None, "metadata": "{}", "tree": None}

###############################################################################

df_entries = pd.read_csv('ABO3_upload_201221.csv')


upload_file_list = ['INCAR','KPOINTS','CONTCAR','OUTCAR','vasprun.xml']

adsorbates = ['OOH', 'O', 'OH']

# exception_list = [12, 15, 25, 31, 34, 59, 67, 68, 77, 78]

###############################################################################
    

with open('api_upload_surf.log', 'w') as f: 
    start_time = time.time()
    
    for idx, index in enumerate(df_entries.index_ori):
        formula = df_entries['formula'][idx]               
        
        try:
            if os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface_np/cont/' % (index, formula)):
                file_path = ('%03d_%s/2nd/surface_np/cont/') % (index, formula)
            elif os.path.exists(os.getcwd() + '/%03d_%s/2nd/surface/2nd/' % (index, formula)):
                file_path = ('%03d_%s/2nd/surface/2nd/') % (index, formula)

            data['tree'] = 'slab/%03d_%s/bare/' % (index, formula)
        
            f.writelines("\nfile_path: %s\n" % file_path)
            f.writelines("upload_path: %s\n" % data['tree'])
        
            for file_name in os.listdir(file_path):
                if file_name in upload_file_list:
                    ret = requests.post(auth['KRICTDB_URL'] + "/api/files/upload",
                                        files = {'file': open(file_path + file_name, "r")}, data = data).text
                    f.writelines([ret,"\n"])

            for k in range(len(adsorbates)):  
                try:
                    if os.path.exists(file_path + '%s/cont/' % adsorbates[k]):
                        if os.path.exists(file_path + '%s/cont/2nd/' % adsorbates[k]):
                            file_path_ads = file_path + '%s/cont/2nd/' % adsorbates[k]
                        else:
                            file_path_ads = file_path + '%s/cont/' % adsorbates[k]
                    elif os.path.exists(file_path + '%s_np/cont/' % adsorbates[k]):
                        if os.path.exists(file_path + '%s_np/cont/2nd/' % adsorbates[k]):
                            file_path_ads = file_path + '%s_np/cont/2nd/' % adsorbates[k]
                        else:
                            file_path_ads = file_path + '%s_np/cont/' % adsorbates[k]
                    elif os.path.exists(file_path + '%s_np/opt/' % adsorbates[k]):
                        file_path_ads = file_path + '%s_np/opt/' % adsorbates[k]
                    else:
                        file_path_ads = file_path + '%s/' % adsorbates[k]
                    print(file_path_ads)
                    

                    data['tree'] = 'slab/%03d_%s/%s/' % (index, formula, adsorbates[k])
                    f.writelines("\tfile_path_ads: %s\n" % file_path_ads)
                    f.writelines("upload_path: %s\n" % data['tree'])

                    for file_name in os.listdir(file_path_ads):
                        if file_name in upload_file_list:
                            ret = requests.post(auth['KRICTDB_URL'] + "/api/files/upload",
                                                files = {'file': open(file_path_ads + file_name, "r")}, data = data).text
                            f.writelines([ret,"\n"])
                except:
                    f.writelines('%03d_%s_%s files are not found !!\n' % (index, formula, adsorbates[k]))

        except:
            f.writelines('%03d_%s files are not found !!\n' % (index, formula))
            continue    

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))
