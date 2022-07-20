# -*- coding: utf-8 -*-
"""
Created on Wed May 18 12:08:23 2022

@author: rbarreror
"""

#
# Import modules
#

import os
import pandas as pd
import numpy as np

from modules.common import explode_be


#
# Test
#

# pdmdf = pd.read_csv(r"C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt", sep="\t")
# pdmdf = pdmdf.loc[:, ['pdm', 'p','q','b','e']].drop_duplicates().reset_index(drop=True)

# og = 'Mouse'

#
# Main
#

def UniProt_gff(pdmdf, og, ph0, qh0, bh0, eh0):
    '''
    Overview:
        Protein sequence annotations are assigned to each peptide. It is
        used gff downloaded from UniProtKB-SwissProt.
        Annotations at peptide level
        
    Input:           
        - pdmdf: Pandas DataFrame
            -p: Peptide sequence
            -q: UniProtKB_AC-ID
            -b: Peptide initial position in q
            -e: Peptide last position in q
            -...
        
        - og: String indicating organism (e.g. Human, Mouse, Pig...)
    
    Output:
        - pdmdf: Pandas DataFrame
            -+=
            -type: Type of annotation
            -start: Start position in protein
            -end: Last position in protein
            -notes: Note assigned to the annotation
            -attr: Additional information
    '''
    
    # Read .feather containing gff annotations
    gff = pd.read_feather(
        os.path.join(
            os.path.dirname(__file__), 
            'data/UniProtKB', 
            og, 
            f'{og}_UniProtKB_SwissProt.feather'
            )
        )
    
    p_ann = pd.merge(
        explode_be(pdmdf).loc[:, ['p','q','b','e']].drop_duplicates(),
        gff,
        how='left',
        left_on='q',
        right_on='UniProtKB_AC-ID'
        )
    
    p_ann = p_ann.loc[
        np.logical_or(
            np.logical_and(
                p_ann['b']>=p_ann['start'],
                p_ann['b']<=p_ann['end']
                ),
            np.logical_and(
                p_ann['e']>=p_ann['start'],
                p_ann['e']<=p_ann['end']
                )
            ),
        :
        ]
    
    p_ann['be_UniProt_gff'] = [
     (int(i), int(j)) 
     for i,j in zip(p_ann['start'].to_list(), p_ann['end'].to_list())
     ]
    
    
    # Add a suffix to extracted column
    colsufix = 'UniProt_gff'
    col = ['type', 'notes', 'attr']
    col = {i:'{i}_{colsufix}' for i in col}
    
    pdmdf = pd.merge(
        pdmdf,
        p_ann.loc[
            :, 
            ['p', 'type', 'notes', 'attr', 'be_UniProt_gff']
            ].rename(col).groupby('p').agg(list),
        how='left',
        on='p'
        )
    
    return pdmdf


# pdmdf_out = UniProt_gff(pdmdf, og)