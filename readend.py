#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import pickle

filename = sys.argv[1]
outpath = sys.argv[2]
with open('./annotation/RibORF_transcript_exon_dict.pickle', 'rb') as f:
    exon_dict = pickle.load(f)

id_list = list(exon_dict.keys())
bamfile = pysam.AlignmentFile(filename, 'rb')

meta3_dict = {}
for item in tqdm(range(len(list(id_list)))):
    name = id_list[item]
    chrom = exon_dict[name].iloc[0,1]
    strand = exon_dict[name].iloc[0,5]
    exon = np.zeros(0)
    if strand == '+':
        for i in range(len(exon_dict[name])):
            exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[i,3], exon_dict[name].iloc[i,4]+1)])
        meta3 = np.zeros((len(exon)+1))
        blocks3 = []
        for read in bamfile.fetch(chrom, exon[0], exon[-1]):
            blocks3.append(read.blocks[-1][-1])
        blocks3 = [x+1 for x in blocks3 if x+1 in exon]
        if len(blocks3) == 0:
            pass
        else:
            metaread3 = []
            for i in range(len(blocks3)):
                metapos3 = np.where(exon == blocks3[i])[0][0]
                metaread3.append(metapos3)
            for i in range(len(metaread3)):
                meta3[metaread3[i]] += 1
        meta3_dict[name] = meta3
    else:
        for i in range(len(exon_dict[name])):
            exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[i,3], exon_dict[name].iloc[i,4]+1)[::-1]])
        meta3 = np.zeros((len(exon)+1))
        blocks3 = []
        for read in bamfile.fetch(chrom, exon[-1], exon[0]):
            blocks3.append(read.blocks[0][0])
        blocks3 = [x+1 for x in blocks3 if x+1 in exon]
        if len(blocks3) == 0:
            pass
        else:
            metaread3 = []
            for i in range(len(blocks3)):
                metapos3 = np.where(exon == blocks3[i])[0][0]
                metaread3.append(metapos3)
            for i in range(len(metaread3)):
                meta3[metaread3[i]] += 1
        meta3_dict[name] = meta3

outname = outpath + '_meta3_RibORF_transcript_dict.pickle'
with open(path2, 'wb') as f:
    pickle.dump(meta3_dict, f)
