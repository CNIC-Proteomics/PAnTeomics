# -*- coding: utf-8 -*-
"""
Created on Wed May 18 13:47:04 2022

@author: rbarreror
"""

#
# Import modules
#

import numpy as np
import gzip
import os
import pandas as pd
import re
import sys

sys.path.append(r'S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Proteomics\GroupTools\Prometeo\src\modules')
from MapID import MapID



#
# Local Functions
#

def read_obo(path):
    #f = open(fr'C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\src\data\GOC\go.obo', 'r')
    f = open(path, 'r')
    
    golist = []
    goelem = []
    for line in f:
        l = line.strip()
        r = re.search(r'\[([\w]+)\]', l)
        if r:
            if len(goelem)>0: 
                golist.append(goelem)
            goelem = []
            
            goelem.append(r.groups()[0])
            continue
        
        r = re.search(r'^(name|id|namespace|def): ([^\n!]+)', l)
        if r:
            goelem.append(r.groups())
            continue
    
    if len(goelem)>0: 
        golist.append(goelem)
    
    f.close()
    
    # Filter Terms
    golist = [sorted(i[1:]) for i in golist if i[0] == 'Term']
    
    # Complete missing fields
    fields = ['name', 'id', 'namespace', 'def']
    golist = [
     i if len(i)==4 else 
     sorted([*i, *[(j, '') for j in fields if j not in list(zip(*i))[0]]] )
     for i in golist
     ]
    
    obodf = pd.DataFrame(
        [list(zip(*i))[1] for i in golist],
        columns=list(zip(*golist[0]))[0]
        )
    
    obodf['def'] = [
     re.search(r'\"([^\"]+)\"', i).groups()[0] if i!='' else '' 
     for i in obodf['def'].to_list()
     ]
    
    return obodf


def read_gaf(path):
    '''
    In this case, we have uniprot id
    For Human and Pig
    '''
    # path = fr'C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\src\data\GOC\{og}\goa_{og.lower()}.gaf.gz'
    
    f = gzip.open(path, 'rb')
    rows = [i.decode('utf-8').strip().split('\t') for i in f if i.decode('utf-8')[0]!='!']
    f.close()
    
    gaf = pd.DataFrame(rows).iloc[:, [1,3,4]].rename({1:'q_id', 3:'GO_qualifier', 4:'GO_id'}, axis=1)
    
    return gaf
    

def add_GO_info(obodf, gaf):
    return pd.merge(
        gaf,
        obodf,
        how='left',
        left_on='GO_id',
        right_on='id',
        ).drop('id', axis=1)
    

#
# Main
#

# Read obo
obodf = read_obo(r'S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Proteomics\GroupTools\Prometeo\src\data\GOC\go.obo')


#
# HUMAN OR PIG
#

og = 'Pig'

path = fr'{og}\goa_{og.lower()}.gaf.gz'

gaf = read_gaf(path)
gaf = gaf.rename({'q_id': 'UniProtKB_AC-ID'}, axis=1)

gaf = pd.merge(
    gaf,
    obodf,
    how='left',
    left_on='GO_id',
    right_on='id',
    ).drop('id', axis=1).groupby('UniProtKB_AC-ID').agg(list).reset_index()

gaf.to_feather(
    f'{og}/{og}_GO.feather'
    )



#
# MOUSE
#
og = 'Mouse'

path = fr'{og}\mgi.gaf.gz'
gaf = read_gaf(path)

# mgi to UniProt 
mgi2up = pd.read_csv(
    os.path.join('Mouse', 'gp2protein.mgi'),
    sep='\t',
    names=['MGI', 'UniProtKB_AC-ID']
    )

mgi2up['MGI'] = [':'.join(i.split(':')[-2:]) for i in mgi2up['MGI'].to_list()]
mgi2up['UniProtKB_AC-ID'] = [i.split(':')[1] for i in mgi2up['UniProtKB_AC-ID'].to_list()]

# Add uniprot id to gaf
gaf = pd.merge(
    gaf,
    mgi2up,
    how='left',
    left_on='q_id',
    right_on='MGI'
    )

gaf = gaf.loc[~pd.isna(gaf['MGI']), :].reset_index(drop=True)

# Add GO information
gaf = add_GO_info(obodf, gaf)

gaf = gaf.drop(['q_id', 'MGI'], axis=1)

gaf = gaf.groupby('UniProtKB_AC-ID').agg(list).reset_index()

gaf.to_feather(
    f'{og}/{og}_GO.feather'
    )



#
# ZEBRAFISH OR RAT
#

# og = 'Zebrafish'
# dbfrom='ZFIN'
# path = r'Zebrafish\zfin.gaf.gz'

og = 'Rat'
dbfrom = 'RGD'
path = r'Rat\rgd.gaf.gz'


gaf = read_gaf(path)

zfin2upi = MapID(gaf['q_id'].drop_duplicates().to_list(), dbfrom, 'UniProtKB-Swiss-Prot')

gaf = pd.merge(
    gaf,
    zfin2upi,
    how='inner',
    left_on='q_id',
    right_on=dbfrom
    )

gaf = add_GO_info(obodf, gaf)

gaf = gaf.rename({'UniProtKB-Swiss-Prot': 'UniProtKB_AC-ID'}, axis=1).drop(['q_id', dbfrom], axis=1)

gaf = gaf.groupby('UniProtKB_AC-ID').agg(list).reset_index()

gaf.to_feather(
    f'{og}/{og}_GO.feather'
    )
