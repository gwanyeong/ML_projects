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

df_entries = pd.read_csv('mpid_list_v2.csv')

# Insert the index of files to upload
num_ini = 0
num_fin = 238

upload_file_list = ['INCAR','KPOINTS','CONTCAR','OUTCAR','vasprun.xml']

###############################################################################

with open('api_upload.log', 'w') as f: 
    start_time = time.time()
    
    for idx in range(num_ini, num_fin):
        formula = df_entries['formula'][idx]
        file_path = '%03d_%s/2nd/DOS/' % (idx + 1.0, formula)
        data['tree'] = 'bulk/%03d_%s/' % (idx + 1.0, formula)
        
        f.writelines("\nfile_path: %s\n" % file_path)
        f.writelines("upload_path: %s\n" % data['tree'])
        
        for file_name in os.listdir(file_path):
            if file_name in upload_file_list:
                ret = requests.post(auth['KRICTDB_URL'] + "/api/files/upload",
                                    files = {'file': open(file_path + file_name, "r")}, data = data).text
                f.writelines([ret,"\n"])

    end_time = time.time()
    f.writelines('Execution time for script (sec) : %6.1f\n' % (end_time - start_time))

