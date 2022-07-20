import gzip
import pandas as pd
import numpy as np

# What is going to be parsed
ga = 'Sscrofa11.1'
filepath = f'{ga}/appris_method.spade.gtf.gz'


# Open gtf.gz
f = gzip.open(filepath, 'rb')

gtf, domType = zip(*[
    [sorted([
        j.strip().split(' ')
        for j in line.decode('utf-8').strip().split('\t')[-1].split(';') # loop subfields in last field
    ]),
    line.decode('utf-8').strip().split('\t')[2]]
    for line in f # loop lines
])

f.close()


gtf = [
    ([j1, j3, *[k.split(':') for k in j2[1].replace('"','').split(',')]])
    for j1,j2,j3 in gtf
]


# Some domains do not have evalue
gtf = [i if 'evalue' in list(zip(*i))[0] else i + [['evalue', np.nan]] for i in gtf]

# Sort fields
gtf = [sorted(i) for i in gtf]
gtf_rows = [list(zip(*i))[1] for i in gtf]
gtf_columns = list(zip(*gtf[0]))[0]


dfgtf = pd.DataFrame(
    [list(zip(*i))[1] for i in gtf], 
    columns=('evalue', 'gene_id', 'hmm_name', 'pep_end', 'pep_start', 'transcript_id')
).astype({'evalue':'float64', 'pep_start':'int64','pep_end':'int64'})

dfgtf['gene_id'] = dfgtf['gene_id'].str.replace('"', '')
dfgtf['transcript_id'] = dfgtf['transcript_id'].str.replace('"', '')

dfgtf['domain_state'] = domType

dfgtf.to_feather(f'{ga}/appris_method.spade.feather')

#dfgtf.to_csv(f'{ga}/appris_method.spade.csv')