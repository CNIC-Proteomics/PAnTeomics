# -*- coding: utf-8 -*-
"""
Created on Fri May 13 12:36:39 2022

@author: rbarreror
"""


#
# IMPORT MODULES
#

import datetime
import pandas as pd
import requests


#
# CONSTANTS
#

dbTable_name = 'dbPTM_table.txt'
download = False


#
# LOCAL FUNCTIONS
#

def download_url(url, save_path, chunk_size=128):
    r = requests.get(url, stream=True)
    with open(save_path, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)

def download_dbPTM(mods):
    _ = [ 
     download_url(
         f'https://awi.cuhk.edu.cn/dbPTM/download/experiment/{i}.zip', 
         f'data/{i}.zip'
         )
     for i in mods
     ]

#
# MAIN
#

if __name__ == '__main__':

    # Read name of the modifications
    f = open('dbPTM_table.txt', 'r')
    mods = [i.strip().split('\t')[0] for i in f]
    f.close()
    
    
    # Download zips
    if download:
        download_dbPTM(mods)
    
    # Read zip tables
    dbptmdf = pd.concat([
     pd.read_csv(f'data/{i}.zip', sep='\t', compression='zip', header=None) for i in mods
    ]).reset_index(drop=True)
    dbptmdf.columns = ['GN', 'UniProtKB_AC-ID', 'n_dbptm', 'mod_dbptm', 'chr_dbptm', 'seq_dbptm']
    
    # Add column containing modified residue
    dbptmdf['modr_dbptm'] = [i[10] if not pd.isna(i) else '' for i in dbptmdf['seq_dbptm'].to_list()]
    
    dbptmdf['chr_dbptm'] = [str(i) for i in dbptmdf['chr_dbptm'].to_list()]
    
    # Save table
    dbptmdf.to_feather('dbPTM.feather')
    dbptmdf.to_csv('dbPTM.tsv', sep='\t', index=False)
    
    # Log
    f = open('dbPTM.log', 'a')
    f.write(f'Update: {datetime.datetime.now().__str__()}\n')
    f.close()
