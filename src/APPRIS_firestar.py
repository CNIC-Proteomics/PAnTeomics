# -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:58:29 2022

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
# Testing variables
#

# pdmdf = pd.read_csv(r"C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt", sep="\t")
# pdmdf = pdmdf.loc[:, ['p','q','b','e']].drop_duplicates()

# from modules.MapID import MapID
# mapdf = MapID(pdmdf['q'].to_list(), 'UniProtKB_AC-ID', 'Ensembl_Transcript')

# ga = 'GRCm39'



#
# Main Function
#

def APPRIS_firestar(pdmdf, mapdf, ga, ph0='p', qh0='q', bh0='b', eh0='e'):
    '''
    Overview:
        Functional residues are assigned to each pdm using APPRIS firestar method
        (cached data).
        Annotations at peptide level
        
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
            -pep_position: Position of the functional residue (in the protein)
            -pep_m: Position in the peptide
            -pep_aa: Residue
            -ligand
    '''
    
    # Read .feather containing functional residues information
    resdf = pd.read_feather(
        os.path.join(
            os.path.dirname(__file__), 
            'data/APPRIS', 
            ga, 
            'appris_method.firestar.feather'
            )
        )
    
    
    # Explode pdm by (b,e) pair
    pdm_be = explode_be(pdmdf)
    
    # Create df with pdm and candidate residues
    pdm_res = pd.merge(
        pdm_be,
        pd.merge(
            mapdf,
            resdf,
            how='left',
            left_on='Ensembl_Transcript',
            right_on='transcript_id'
            ),
        left_on='q',
        right_on='UniProtKB_AC-ID'
        )
    
    # Create df with pdm and filtered residue
    pdm_res = pdm_res.loc[
        np.logical_and(
            pdm_res['pep_position'] > pdm_res['b'],
            pdm_res['pep_position'] < pdm_res['e']    
            ),
        :
        ].reset_index(drop=True)
    
        
    # Get position (in p) and aa of functional residue
    p_res = pdm_res.loc[:, ['p','b','pep_position', 'ligands']].drop_duplicates()
    p_res['pep_m'] = p_res['pep_position']-p_res['b']+1
    p_res['pep_aa'] = [
        i[int(j)] for i,j in zip(p_res['p'].to_list(), p_res['pep_m']-1)
        ]
    
    pdm_res = pd.merge(
        pdm_res,
        p_res.loc[:, ['p', 'b', 'pep_position', 'ligands', 'pep_m', 'pep_aa']],
        how='left',
        on=['p','b','pep_position','ligands']
        )
    
    
    # Add information to pdmdf
    pdmdf_out = pd.merge(
        pdmdf,
        pdm_res.loc[
            :, 
            ['p', 'pep_position', 'ligands', 'pep_m', 'pep_aa']
            ].drop_duplicates().groupby('p').agg(list),
        
        how='left',
        on='p'
        )
    
    # Convert float to int
    pdmdf_out['pep_position'] = [
        [int(j) for j in i] if type(i)==list else i 
        for i in pdmdf_out['pep_position'].to_list()
        ]
    
    pdmdf_out['pep_m'] = [
        [int(j) for j in i] if type(i)==list else i 
        for i in pdmdf_out['pep_m'].to_list()
        ]
    
    return pdmdf_out


# pdmdf_out = APPRIS_firestar(pdmdf, mapdf, ga)