# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:57:15 2022

@author: rbarreror
"""

#
# Import modules
# 

import numpy as np
import os
import pandas as pd

from modules.common import explode_be

#
# Testing variables
#

# pdmdf = pd.read_csv(r"S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Proteomics\GroupTools\PAnTeomics\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt", sep="\t")
# pdmdf = pdmdf.loc[:, ['pdm','p','q','b','e','n', 'm']].drop_duplicates()


#
# Main Function
#

def dbPTM(pdmdf, pdmh0, ph0, qh0, bh0, eh0, mh0):
    '''
    Overview:
        pdm with "modifiable" residues contained in dbPTM are annotated (cached data).
        Annotations at peptide level.
        
    Input:           
        - pdmdf: Pandas DataFrame
            -p: Peptide sequence
            -q: UniProtKB_AC-ID
            -b: Peptide initial position in q
            -e: Peptide last position in q
            -...
        
    Output:
        - pdmdf: Pandas DataFrame
            -+=
            -n_dbptm: Position of the annotated residue in the protein
            -mod_dbptm: Name of the modification
            -modr_dbptm: modifiable residue
            -m_dbptm: Position of the annotated residue in the peptide
    '''
        
    # read dbPTM data
    dbptmdf = pd.read_feather(
        os.path.join(
            #r'C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\src',
            os.path.dirname(__file__),
            'data/dbPTM/dbPTM.feather'
            )
        )
    
    # pdm2p = pdmdf.loc[:, [pdmh0, ph0]]
    
    p_db = pd.merge(
        explode_be(pdmdf.drop([pdmh0, mh0], axis=1).drop_duplicates()),
        dbptmdf,
        how='left',
        left_on=qh0,
        right_on='UniProtKB_AC-ID'
        )
    
    p_db = p_db.loc[
        np.logical_and(
            p_db[bh0]<=p_db['n_dbptm'],
            p_db[eh0]>=p_db['n_dbptm']
            ),
        [ph0, bh0, *dbptmdf.columns.to_list()]
        ]
    
    p_db['m_dbptm'] = (p_db['n_dbptm']-p_db[bh0]+1).astype(int)
    
    p_db = p_db.drop(
        [bh0, 'GN', 'UniProtKB_AC-ID', 'seq_dbptm', 'chr_dbptm'], axis=1
        ).groupby(ph0).agg(list).reset_index()
    
    p_db['n_dbptm'] = [
        [int(j) for j in i] if type(i)==list else i 
        for i in p_db['n_dbptm'].to_list()
        ]
    
    
    pdmdf_out = pd.merge(
        pdmdf,
        p_db,
        how='left',
        on=ph0
        )
    
    pdmdf_out['mod_dbptm_in_m'] = [
     list(set([j3 for j2, j3 in zip(i2, i3) if i1==j2])) if type(i2)==list
     else []
     #i3
     for i1, i2, i3 in zip(
             pdmdf_out[mh0].to_list(),
             pdmdf_out['m_dbptm'].to_list(),
             pdmdf_out['mod_dbptm'].to_list()
             )
    ]
    
    return pdmdf_out

# pdmdf_out = dbPTM(pdmdf, 'pdm', 'p', 'q', 'b', 'e', 'm')
