#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle

os.makedirs('./annotation', exist_ok=True)

gtf = pd.read_csv('./gencode.v31.annotation.gtf', sep='\t', header=None)
gtf_exon = gtf[gtf[2] == 'exon']
ids = [x.split(';')[1].split('"')[1] for x in gtf_exon[8]]
names = [x.split(';')[3].split('"')[1] for x in gtf_exon[8]]
symbol = pd.Series(names, index=ids).drop_duplicates()

res = pd.read_csv('./RibORF_output/pred.pvalue.parameters.txt', sep='\t', index_col=0)
res = res[res['pred.pvalue']>0.7]
res.to_csv('./RibORF_output/RibORF_result.csv')
orfID = res['orfID']
df = pd.concat([orfID.iloc[:,:1], orfID[2].str.split('|', expand=True), \
orfID.iloc[:,3], orfID[4].str.split('|', expand=True)], axis=1)
df.columns = ['tid','chrom','strand','ORF number','tlen','tcoord','tstop','orftype','codon']
df['tfam'] = list(symbol.reindex(df['tid']))
df['tcoord'] = df['tcoord']-1
df['tstop'] = df['tstop']-1
df['frame'] = df['tid']+'_'+(df['tcoord']%3).astype(str)
df.to_pickle('./annotation/orfID.pickle')

transcripts = list(set(df['tid']))
gtf_exon_RibORF = gtf_exon[gtf_exon['transcript_id'].isin(transcripts)]
exon_dict = {}
for name, group in gtf_exon_RibORF.groupby('transcript_id'):
    exon_dict[name] = group
    exon_dict[name] = exon_dict[name].reset_index()
with open('./annotation/RibORF_transcript_exon_dict.pickle', 'wb') as f:
    pickle.dump(exon_dict, f)

fasta_in = './gencode.v31.transcripts.fa'
seq_dict = {}
for record in SeqIO.parse(fasta_in, 'fasta'):
    id_part = record.id
    name = id_part.split('|')[0]
    seq = record.seq
    if name in transcripts:
        seq_dict[id_part] = seq
with open('./annotation/RibORF_transcripts_sequence_dict.pickle', 'wb') as f:
    pickle.dump(seq_dict, f)
