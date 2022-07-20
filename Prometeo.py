# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:05:17 2022

@author: rbarreror
"""

import sys
sys.path.insert(0,r'C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\src')


import pandas as pd

from src.APPRIS_firestar import APPRIS_firestar
from src.APPRIS_spade import APPRIS_spade

#
# Testing variables
#

pdmdf = pd.read_csv(r"C:\Users\rbarreror\home\Proteomics\GroupTools\Prometeo\test\Heteroplasmy\PDMTable_1_Heart_PDMTable1.txt", sep="\t")
pdmdf = pdmdf.loc[:, ['pdm','p','q','b','e']]

from src.modules.MapID import MapID
mapdf = MapID(pdmdf['q'].to_list(), 'UniProtKB_AC-ID', 'Ensembl_Transcript')

ga = 'GRCm39'


APPRIS_firestar(pdmdf, mapdf, ga)
APPRIS_spade(pdmdf, mapdf, ga)
