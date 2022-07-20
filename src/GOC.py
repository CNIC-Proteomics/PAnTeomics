# -*- coding: utf-8 -*-
"""
Created on Fri May 20 09:44:32 2022

@author: rbarreror
"""


#
# Import Modules
#

import os
import pandas as pd

# 
# Test
#

# pdmdf = pd.read_csv(r"S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Proteomics\GroupTools\Prometeo\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt", sep="\t")
# pdmdf = pdmdf.loc[:, ['pdm', 'p','q','b','e']].drop_duplicates().reset_index(drop=True)

# og = 'Mouse'


#
# Main
#

def GOC(pdmdf, og, qh0):
    '''
    Overview:
        GO annotations are assigned to each protein. 
        Obtained from .GAF files downloaded from:
            http://current.geneontology.org/products/pages/downloads.html
        Annotations at protein level
        
    Input:           
        - pdmdf: Pandas DataFrame
            -q: UniProtKB_AC-ID
            -og: String indicating the organism (e.g. Human, Mouse...)
        
        - og: String indicating organism (e.g. Human, Mouse, Pig...)
    
    Output:
        - pdmdf: Pandas DataFrame
            -+=
            -GO_qualifier: Relation between the protein and the GO term. E.g. involved in
            -GO_id: E.g. GO:0003677
            -def: E.g. Cytoplasm situated near, or occurring around, the nucleus.
            -name: E.g. perinuclear region of cytoplasm
            -namespace: cellular_component
    '''
    

    # Read feather containing GO information
    gocdf = pd.read_feather(
        os.path.join(
            os.path.dirname(__file__),
            'data/GOC',
            og,
            f'{og}_GO.feather'
            )
        )
    
    pdmdf_out = pd.merge(
        pdmdf,
        gocdf.rename({'UniProtKB_AC-ID':qh0}, axis=1),
        how='left',
        on=qh0,
        )
    
    return pdmdf_out


#pdmdf_out = GOC(pdmdf, og)