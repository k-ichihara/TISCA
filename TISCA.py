#!/usr/bin/env python
import pandas as pd
import numpy as np
import pickle
import pysam
from tqdm import tqdm
from Bio import SeqIO
import os
from scipy.signal import argrelmax
import csv

#Read end TIS peak and EIF3 scanning release
scan_left = 22
scan_right = 26
#Scanning release peak order
order1 = 1
width = 15
gap = 10

C_in = 'Ribo'
L_in = 'GTI'
S_in = 'Sel'

orfID = pd.read_pickle('./data/orfID.pickle')
Ribo = orfID.set_index('tid')['tcoord']

with open('./annotation/RibORF_transcripts_sequence_dict.pickle', 'rb') as f:
    seq_dict = pickle.load(f)
with open('./annotation/RibORF_transcript_dict.pickle', 'rb') as f:
    exon_dict = pickle.load(f)
id_list = list(exon_dict.keys())

def bam_count(inname):
    bamfile = pysam.AlignmentFile(inname, "rb")
    return bamfile.count()

def in_frame(seq):
    if len(seq) % 3 ==1:
        seq += 'NN'
    elif len(seq) % 3 == 2:
        seq += 'N'
    else:
        pass
    return seq

def add_position(df):
    pos_list = []
    for i in range(len(df)):
        name = df['tid'][i]
        TIS_start = df['TIS peak start'][i] + df['TIS offset'][i]
        chrom = exon_dict[name].iloc[0,1]
        strand = exon_dict[name].iloc[0,5]
        exon = np.zeros(0)
        if strand == '+':
            for j in range(len(exon_dict[name])):
                exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[j,3], exon_dict[name].iloc[j,4]+1)])
            pos_list.append(chrom + ':' + str(int(exon[TIS_start])) + '-' + str(int(exon[TIS_start+1])))
        else:
            for j in range(len(exon_dict[name])):
                exon = np.concatenate([exon, np.arange(exon_dict[name].iloc[j,3], exon_dict[name].iloc[j,4]+1)[::-1]])
            pos_list.append(chrom + ':' + str(int(exon[TIS_start])) + '-' + str(int(exon[TIS_start-1])))
    df['Position'] = pos_list

def scanning_release(Pro_in, order, width, gap):
    with open('./annotation/' + Pro_in + '/' + Pro_in + '_RibORF_transcript_meta_dict.pickle', 'rb') as f:
        read_dict = pickle.load(f)
    with open('./annotation/' + Pro_in + '/' + Pro_in + '_meta3_RibORF_transcript_dict.pickle', 'rb') as f:
        end_dict = pickle.load(f)
    meta3_peak_dict = {}
    for item in range(len(list(id_list))):
        name = id_list[item]
        y = end_dict[name]
        peaks = np.stack([argrelmax(y, order = order)[0], y[argrelmax(y, order = order)]])
        peaks = peaks[:, peaks[1,:] > 0]
        if peaks.sum() == 0:
            pass
        else:
            for i, item in enumerate(peaks[0,:]):
                left_start = max([0, item-gap-width])
                left_stop = max([0, item-gap])
                right_start = min([len(read_dict[name]), item+gap])
                right_stop = min([len(read_dict[name]), item+gap+width])
                left = read_dict[name][int(left_start):int(left_stop)]
                left = left[left >= 3]
                right = read_dict[name][int(right_start):int(right_stop)]
                peaks[1,i] = max([1, left.sum()]) / max([1, right.sum()])
            peaks = peaks[:, peaks[1,:] > 1]
            if peaks.sum() == 0:
                pass
            else:
                meta3_peak_dict[name] = peaks
    return meta3_peak_dict

def TIS_peak(LTM_in):
    peak_region_dict = {}
    with open('./annotation/' + LTM_in + '/' + LTM_in + '_RibORF_transcript_meta_dict.pickle', 'rb') as f:
        LTM_meta_transcript_dict = pickle.load(f)
        LTM_tot = bam_count("./data/" + LTM_in + ".bam")
    for item in range(len(id_list)):
        name = id_list[item]
        LTM = LTM_meta_transcript_dict[name] * 1000000 / LTM_tot
        LTM_diff = np.gradient(LTM)
        std = np.std(LTM_diff)
        out = np.where(LTM_diff >= std*4)[0]

        region_name = []
        x_plus = out
        if len(out) <= 1:
            pass
        else:
            div = np.where(np.diff(out) > 1)[0] + 1
            if len(div) == 1:
                if len(out[:div[0]]) > 1:
                    region_name.append(out[:div[0]])
            elif len(div) > 1:
                if len(out[:div[0]]) > 1:
                    region_name.append(out[:div[0]])
                for i in range(1, len(div)):
                    region_name_div = out[div[i-1]:div[i]]
                    if len(region_name_div) < 2:
                        pass
                    else:
                        region_name.append(region_name_div)
                if len(x_plus[div[-1]:]) > 1:
                    region_name.append(out[div[-1]:])
            else:
                region_name.append(out)
        peaks = [x[1] for x in region_name]
        if len(peaks) > 0:
            peak_region_dict[name] = [[x, LTM[x:x+30].sum()] for x in peaks]

    return peak_region_dict

def TIS_scanning(CHX_in, LTM_in, Pro_in):
    id_list = list(set(release_dict.keys()) & set(TIS_dict.keys()))
    res = ['tid', 'TIS peak start', 'TIS peak area', 'Scanning release point', 'Scanning release score',
            'TIS offset', 'TIS start codon', 'TIS ORF sequence']
    peak3_eff_ORF_dict = {}
    for item in tqdm(range(len(id_list))):
        name = id_list[item]
        peak3_eff_ORF_name = []
        peak_start = [x[0] for x in TIS_dict[name]]
        peak_area = [x[1] for x in TIS_dict[name]]
        scan_point = release_dict[name]
        Ribo_name = Ribo[Ribo.index == name]
        for j in range(len(peak_start)):
            for offset in range(11,14):
                TIS_start_offset = int(peak_start[j] + offset)
                peak3_eff = scan_point[:, (scan_point[0] <= TIS_start_offset + scan_right) & (scan_point[0] >= TIS_start_offset + scan_left)]
                if peak3_eff.sum() == 0:
                    pass
                else:
                    start_codon = str(seq_dict[name][TIS_start_offset:TIS_start_offset+3])
                    for i in range(len(Ribo_name)):
                        ORF_start = Ribo_name.values[i]
                        if ORF_start == TIS_start_offset:
                            release_point = int(peak3_eff[0][0])
                            release_score = peak3_eff[1][0]
                            start_codon_ex1 = str(seq_dict[name][TIS_start_offset-1:TIS_start_offset+4])
                            start_codon_ex2 = str(seq_dict[name][TIS_start_offset-2:TIS_start_offset+5])
                            TIS_seq = seq_dict[name][TIS_start_offset:]
                            TIS_ORF = in_frame(TIS_seq).translate(to_stop=True) + '*'
                            peak3_eff_ORF_name.append([name, peak_start[j], peak_area[j], release_point, release_score, offset,
                                                       start_codon, str(TIS_ORF)])
        if not peak3_eff_ORF_name:
            pass
        else:
            peak3_eff_ORF_df_name = pd.DataFrame(peak3_eff_ORF_name)
            peak3_eff_ORF_df_name.columns = res
            peak3_eff_ORF_dict[name] = peak3_eff_ORF_df_name

    peak3_eff_ORF_df = pd.concat([v for v in peak3_eff_ORF_dict.values()])
    peak3_eff_ORF_df = peak3_eff_ORF_df.reset_index().iloc[:,1:]
    peak3_eff_ORF_df.index = peak3_eff_ORF_df['tid'] + '_' + peak3_eff_ORF_df['TIS peak start'].astype(str) + '_' + peak3_eff_ORF_df['TIS offset'].astype(str)
    add_position(peak3_eff_ORF_df)
    peak3_eff_ORF_df['TIS_ID'] = peak3_eff_ORF_df['tid'] + '_' + peak3_eff_ORF_df['TIS peak start'].astype(str)
    return peak3_eff_ORF_df

# Main
filename_GTI = L_in + '-' + C_in
TIS_dict = TIS_peak(L_in)
release_dict = scanning_release(S_in, order1, width, gap)
df = TIS_scanning(C_in, L_in, S_in)
df.to_csv('./TISCA_result_'+ filename_GTI + '_' + str(scan_left) + '-' + str(scan_right) + '_' + S_in + '.csv')
