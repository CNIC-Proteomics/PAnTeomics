# -*- coding: utf-8 -*-
"""
Created on Fri May 27 13:01:33 2022

@author: rbarreror
"""

import pandas as pd

from modules.MapID import MapID
from APPRIS_spade import APPRIS_spade
from APPRIS_firestar import APPRIS_firestar
from dbPTM import dbPTM
from UniProt_gff import UniProt_gff
from GOC import GOC



#repdf = pd.read_csv(r"S:\U_Proteomica\UNIDAD\DatosCrudos\Jose_Antonio_Enriquez\Ratones_Heteroplasmicos\ReprocesadoRafa\BorradorPaper\Procesamiento_CD\iSanXoT\Heart-1\reports\pdm_table.tsv", sep="\t", header=[0,1], low_memory=False)
# #pdmdf = pdmdf.loc[:, ['pdm', 'p','q','b','e']].drop_duplicates().reset_index(drop=True)

#
# Variables
#

infile = r"S:\U_Proteomica\LABS\LAB_ARR\LaCaixa\tejidos-secretomas\New_Comet_500\iSanxot_Intima\QUANTWF2\reports\NM\LIMMA_NM_pgmpq_table.tsv"
outfile = r'S:\U_Proteomica\LABS\LAB_ARR\LaCaixa\tejidos-secretomas\New_Comet_500\iSanxot_Intima\QUANTWF2\reports\NM\LIMMA_NM_pgmpq_table_annotated2.tsv'

og = 'Human'

pdmh = ('pdm', 'LEVEL')
ph = ('p', 'REL')
qh = ('q', 'REL')
bh = ('b', 'REL')
eh = ('e', 'REL')
mh = ('m', 'REL')

mapdf = None

#
# End Variables
#





def PAnTeomics(repdf, pdmh, ph, qh, bh, eh, mh, og, mapdf=None):

    #
    # Constants
    #
    
    idx = pd.IndexSlice
    
    # Appris works with genome assembly
    og2ga ={
            'Human':'GRCh38',
            'Mouse':'GRCm39',
            'Pig':'Sscrofa11.1',
            'Rat':'mRatBN7.2',
            'Zebrafish':'GRCz11',
            }
    
    
    if not mapdf:
        mapdf = MapID(repdf[qh[0], qh[1]].drop_duplicates().tolist(), 'UniProtKB_AC-ID', 'Ensembl_Transcript')
    
    pdmdf_pdm_input = repdf.loc[
        :,
        idx[[pdmh[0],ph[0],qh[0],bh[0],eh[0],mh[0]], :]
        ].droplevel(1, axis=1).astype({bh[0]:str, eh[0]:str}).drop_duplicates()
    
    pdmdf_p_input = repdf.loc[
        :,
        idx[[ph[0],qh[0],bh[0],eh[0]], :]
        ].droplevel(1, axis=1).astype({bh[0]:str, eh[0]:str}).drop_duplicates()
    
    pdmdf_q_input = repdf.loc[
        :,
        idx[[qh[0]], :]
        ].droplevel(1, axis=1).drop_duplicates()
    
    
    #
    # APPRIS SPADE
    #
    
    tmp = APPRIS_spade(
        pdmdf_p_input,
        mapdf,
        og2ga[og],
        ph[0], qh[0], bh[0], eh[0]
    ).loc[:, [ph[0], qh[0], 'hmm_name', 'be_domain', 'domain_state']]
    
    tmp.columns = pd.MultiIndex.from_product([tmp.columns, ['APPRIS_SPADE']])
    
    repdf = pd.merge(
        repdf,
        tmp,
        left_on=[ph, qh],
        right_on=[(ph[0], 'APPRIS_SPADE'), (qh[0], 'APPRIS_SPADE')],
        how='left'
    ).drop([(ph[0], 'APPRIS_SPADE'), (qh[0], 'APPRIS_SPADE')], axis=1)
    
    
    
    #
    # APPRIS FIRESTAR
    #
    
    tmp = APPRIS_firestar(
        pdmdf_p_input,
        mapdf,
        og2ga[og],
        ph[0], qh[0], bh[0], eh[0]
    ).loc[:, [ph[0], qh[0], 'pep_position', 'ligands', 'pep_m', 'pep_aa']]
    
    tmp.columns = pd.MultiIndex.from_product([tmp.columns, ['APPRIS_FIRESTAR']])
    
    repdf = pd.merge(
        repdf,
        tmp,
        left_on=[ph, qh],
        right_on=[(ph[0], 'APPRIS_FIRESTAR'), (qh[0], 'APPRIS_FIRESTAR')],
        how='left'
    ).drop([(ph[0], 'APPRIS_FIRESTAR'), (qh[0], 'APPRIS_FIRESTAR')], axis=1)
    
    
    #
    # dbPTM
    #
    
    tmp = dbPTM(
        pdmdf_pdm_input, pdmh[0], ph[0], qh[0], bh[0], eh[0], mh[0]
    ).loc[:, [pdmh[0], 'mod_dbptm_in_m']]# 'n_dbptm', 'mod_dbptm', 'modr_dbptm', 'm_dbptm']]
    
    tmp.columns = pd.MultiIndex.from_product([tmp.columns, ['dbPTM']])
    
    repdf = pd.merge(
        repdf,
        tmp,
        left_on=[pdmh],
        right_on=[(pdmh[0], 'dbPTM')],
        how='left'
    ).drop([(pdmh[0], 'dbPTM')], axis=1)
    
    
    
        
    #
    # UniProt_gff
    #
    
    tmp = UniProt_gff(
        pdmdf_p_input,
        og,
        ph[0], qh[0], bh[0], eh[0]
    ).loc[:, [ph[0],qh[0],'type', 'notes', 'attr', 'be_UniProt_gff']]
    
    tmp.columns = pd.MultiIndex.from_product([tmp.columns, ['UniProtKB']])
    
    repdf = pd.merge(
        repdf,
        tmp,
        left_on=[ph, qh],
        right_on=[(ph[0], 'UniProtKB'), (qh[0], 'UniProtKB')],
        how='left'
    ).drop([(ph[0], 'UniProtKB'), (qh[0], 'UniProtKB')], axis=1)
    
    
    #
    # GOC
    #
    
    tmp = GOC(
        pdmdf_q_input,
        og,
        qh[0]
    ).loc[:, [qh[0],'GO_qualifier', 'GO_id', 'def', 'name','namespace']]
    
    tmp.columns = pd.MultiIndex.from_product([tmp.columns, ['GOC']])
    
    repdf = pd.merge(
        repdf,
        tmp,
        left_on=[(qh[0], 'REL')],
        right_on=[(qh[0], 'GOC')],
        how='left'
    ).drop([(qh[0], 'GOC')], axis=1)
    
    
    return repdf


# if __name__ == '__main__':

# Read report
repdf = pd.read_csv(infile, sep="\t", header=[0,1], low_memory=False)

# annotate
repdf_out = PAnTeomics(repdf, pdmh, ph, qh, bh, eh, mh, og, mapdf)

# Write output
repdf_out.drop(
    ['GO_qualifier', 'def', 'name', 'namespace'], axis=1, level=0
    ).to_csv(outfile, sep='\t', index=False)

#repdf_out.to_csv(r'S:\U_Proteomica\LABS\LAB_ARR\LaCaixa\tejidos-secretomas\New_Comet_500\iSanxot_Intima\QUANTWF2\reports\NM\LIMMA_NM_pgmpq_table_annotated.tsv')
