# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:39:19 2022

@author: rbarreror
"""

#
# Import Modules
#

import numpy as np
import os
import pandas as pd

from modules.common import explode_be


#
# Testing module
# 


# pdmTable_path = r'C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt'
# pdmTable = pd.read_csv(pdmTable_path, sep='\t')
# pdmdf = pdmTable.loc[:, ['p', 'q', 'b', 'e']].drop_duplicates()

# # UniProt -> Ensemble_transcript
# from modules.MapID import MapID
# mapdf = MapID(pdmdf['q'], 'UniProtKB_AC-ID', 'Ensembl_Transcript')

# # Gene Assembly
# ga = 'GRCm39'



#
# Main Function
#

def APPRIS_spade(pdmdf, mapdf, ga, ph0='p', qh0='q', bh0='b', eh0='e'):
    '''
    Overview:
        Protein domains are assigned to each p using APPRIS spade method
        (cached data).
        Annotations at peptide level.
        
    Input:           
        - pdmdf: Pandas DataFrame
            -p: Peptide sequence
            -q: UniProtKB_AC-ID
            -b: Peptide initial position in q
            -e: Peptide last position in q
            -...
        
        - mapdf: Pandas DataFrame
            -from: UniProtKB_AC-ID
            -to: Ensembl_Transcript
        
        - ga: String indicating Genome Assembly (e.g. GRCh38, GRCm39...)
    
    Output:
        - pdmdf: Pandas DataFrame
            -+=
            -hmm_name: Domain name
            -be_domain: Begin and End of domain in protein
    '''
    
    # Read dataframe with Ensemble_Transcript -> Domain
    domdf = pd.read_feather(
        os.path.join(
            os.path.dirname(__file__),
            'data/APPRIS', 
            ga, 
            'appris_method.spade.feather'
            )
        )
    
    # UniProt -> Domain
    q_dom_df = pd.merge(
        mapdf,
        domdf,
        how='left',
        left_on='Ensembl_Transcript',
        right_on='transcript_id'
    )
    
    # pdm,UniProt,Pos
    pdm_be = explode_be(pdmdf)
    
    
    # pdm,UniProt,Pos,Domain 
    pdm_qpos_dom = pd.merge(
        pdm_be,
        q_dom_df,
        how='left',
        left_on=qh0,
        right_on='UniProtKB_AC-ID'
    )
    
    
    # check if b or e is in (pep_start, pep_end) interval (it means intersection)
    pdm_qpos_dom['intersect'] = np.logical_or(
        np.logical_and(
            pdm_qpos_dom[bh0]>pdm_qpos_dom['pep_start'],
            pdm_qpos_dom[bh0]<pdm_qpos_dom['pep_end']
        ),
        np.logical_and(
            pdm_qpos_dom[eh0]>pdm_qpos_dom['pep_start'],
            pdm_qpos_dom[eh0]<pdm_qpos_dom['pep_end']
        )
    )
    
    # A column with (b,e) of domain
    pdm_qpos_dom['be_domain'] = [
        (i,j) 
        for i,j in zip(
                pdm_qpos_dom['pep_start'].to_list(), 
                pdm_qpos_dom['pep_end'].to_list()
                )
        ]
    
    # Filter pdm,UniProt,Pos,Domain with intersect
    pdm_dom = pdm_qpos_dom.loc[
        pdm_qpos_dom['intersect'], 
        [ph0, 'hmm_name', 'be_domain', 'domain_state']
        ].drop_duplicates().groupby([ph0]).agg(list).reset_index()
    
    pdm_dom['be_domain'] = [
        [tuple([int(j[0]), int(j[1])]) for j in i] if type(i)==list else i 
        for i in pdm_dom['be_domain'].to_list()
        ]
    
    
    pdmdf_out = pd.merge(
        pdmdf,
        pdm_dom,
        how='left',
        on=ph0
    )

    return pdmdf_out


# pdmdf_out = APPRIS_spade(pdmdf, mapdf, ga)