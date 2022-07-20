# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:20:21 2022

@author: rbarreror
"""


#
# Import modules
#

import gzip
import os
import pandas as pd


#
# Set constants
#

og = 'Zebrafish'


#
# Main
#

path = os.path.join(og, f'{og}_UniProtKB_SwissProt.gff.gz')

# Get readable lines
f = gzip.open(path, 'rb')
lines = [i.decode('utf-8').strip().split('\t') for i in f]
lines = [i for i in lines if len(i)==9]
f.close()

# build pandas dataframe
gff = pd.DataFrame(
    lines,
    columns=[
        'UniProtKB_AC-ID',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'attr'
        ]
    )

# Extract notes
notes = [
 [
  j.split('=')[1] 
  for j in i.split(';') if j.split('=')[0] == 'Note'
  ] 
 for i in gff['attr'].to_list()
 ]

notes = ['' if len(i)==0 else i[0] for i in notes]

gff['notes'] = notes

gff = gff.astype({'start': 'int64', 'end': 'int64'})

# Write feather
gff.loc[
    :,
    ['UniProtKB_AC-ID', 'type', 'start', 'end', 'notes', 'attr']
    ].to_feather(
    path.split('.')[0]+'.feather'
    )


