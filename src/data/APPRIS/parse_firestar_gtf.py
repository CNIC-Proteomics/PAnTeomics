# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:16:05 2022

@author: rbarreror
"""

import pandas as pd

ga = 'Sscrofa11.1'

filepath = f'{ga}/appris_method.firestar.gtf.gz'

df = pd.read_csv(filepath, sep='\t', header=None, compression='gzip', low_memory=False)

info = [
 sorted(
        [
            j.strip().replace('"', '').split(' ') # split subfield
            for j in i.split(';') # loop subfields
        ]
        ) 
 for i in df.iloc[:, -1].to_list() # loop elements in last column
]

info = [[i1,i3,*[k.split(':') for k in i2[1].split(',')]] for i1,i2,i3 in info]

columns = list(zip(*info[0]))[0]
rows = [list(zip(*i))[1] for i in info]

infodf = pd.DataFrame(
    rows,
    columns=columns
    ).astype({'pep_position': 'int64'})

# write output
infodf.to_feather(f'{ga}/appris_method.firestar.feather')
