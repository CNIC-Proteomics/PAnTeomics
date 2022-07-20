# -*- coding: utf-8 -*-
"""
Created on Wed May 11 10:46:46 2022

@author: rbarreror
"""

def explode_be(pdmdf):
    '''
    Overview:
        A pdm can appear many times in the protein, yielding multiple b:e
        pairs. This function explode a pdm by b:e set.
        Input dataframe must contain b and e column, with values
        separated by ;
    '''
    #pdm_qpos = pdmdf.loc[:, ['pdm', 'p', 'q']]
    pdm_qpos = pdmdf.copy()
    pdm_qpos['be'] = [list(zip(*[i.split(';'),j.split(';')])) for i,j in zip(pdmdf['b'].to_list(), pdmdf['e'].to_list())]
    pdm_qpos = pdm_qpos.explode('be').reset_index(drop=True)
    pdm_qpos['b'] = list(zip(*pdm_qpos['be'].to_list()))[0]
    pdm_qpos['e'] = list(zip(*pdm_qpos['be'].to_list()))[1]
    pdm_qpos = pdm_qpos.astype({'b':'int64', 'e':'int64'})
    
    return pdm_qpos.loc[:, pdmdf.columns]


def batch(iterable, n=1):
    '''
    Split iterable in chunks of equal length
    https://stackoverflow.com/questions/8290397/how-to-split-an-iterable-in-constant-size-chunks
    '''
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]