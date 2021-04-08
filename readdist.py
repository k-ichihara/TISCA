#!/usr/bin/env python
import os
import cys
import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import pickle
import itertools

filename = sys.argv[1]
outpath = sys.argv[2]
with open('./annotation/RibORF_transcript_exon_dict.pickle', 'rb') as f:
    exon_dict = pickle.load(f)

id_list = list(exon_dict.keys())
bamfile = pysam.AlignmentFile(filename, 'rb')

meta1_dict = {}
for item in tqdm(range(len(list(id_list)))):
    name = id_list[item]
    chrom = exon_dict[name].iloc[0,1]
    strand = exon_dict[name].iloc[0,5]
    exon = np.zeros(0)
    if strand == '+':
        for i in range(len(exon_dict[name])):
            exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[i,3], exon_dict[name].iloc[i,4]+1)])
        meta1 = np.zeros((len(exon)+1))
        read_position = []
        for read in bamfile.fetch(chrom, exon[0], exon[-1]):
            read_position.append(read.get_reference_positions())
        read_position_flat = list(itertools.chain.from_iterable(read_position))
        read_position_flat = [x+1 for x in read_position_flat if x+1 in exon]
        metaread = []
        for i in range(len(read_position_flat)):
            metapos = np.where(exon == read_position_flat[i])[0][0]
            metaread.append(metapos)
        for i in range(len(metaread)):
            meta1[metaread[i]] += 1
        meta1_dict[name] = meta1
    else:
        for i in range(len(exon_dict[name])):
            exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[i,3], exon_dict[name].iloc[i,4]+1)[::-1]])
        meta1 = np.zeros((len(exon)+1))
        read_position = []
        for read in bamfile.fetch(chrom, exon[-1], exon[0]):
            read_position.append(read.get_reference_positions())
        read_position_flat = list(itertools.chain.from_iterable(read_position))
        read_position_flat = [x+1 for x in read_position_flat if x+1 in exon]
        metaread = []
        for i in range(len(read_position_flat)):
            metapos = np.where(exon == read_position_flat[i])[0][0]
            metaread.append(metapos)
        for i in range(len(metaread)):
            meta1[metaread[i]] += 1
        meta1_dict[name] = meta1

outname = outpath + '_RibORF_transcript_meta_dict.pickle'
with open(outpath, 'wb') as f:
    pickle.dump(meta1_dict, f)
