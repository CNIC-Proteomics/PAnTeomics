# -*- coding: utf-8 -*-
"""
Created on Wed May 11 10:25:15 2022

@author: rbarreror
"""


#
# Import Modules
#


import pandas as pd
import requests
import time


#
# Functions
# 

def MapID(ids, fromdb, todb):
    '''
    '''
    # Send request
    ids = ','.join(ids)

    url_map = 'https://rest.uniprot.org/idmapping/run'

    r_url = requests.post(url_map, {
        'from': fromdb,
        'to': todb,
        'ids': ids,
    })

    if not r_url.ok:
        return pd.DataFrame(columns=[fromdb,fromdb])
        
        
    # Get Ensembl id
    while True:
        url_get = f'https://rest.uniprot.org/idmapping/stream/{r_url.json()["jobId"]}?compressed=false&format=json'
        r_mapped = requests.get(url_get)
        
        if r_mapped.ok:
            r_mapped = r_mapped.json()['results']
            break
        else:
            #print('** Waiting map...')
            time.sleep(5)

    # Create dataframe mapping UniProt to Ensemble
    mapdf = pd.DataFrame(
        [(i['from'], i['to']) for i in r_mapped],
        columns=[fromdb, todb]
    )
    
    return mapdf
        