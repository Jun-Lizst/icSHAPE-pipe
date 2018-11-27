#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import re
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import pandas as pd
import version

def gene_type(raw_type):
    valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA', 'mRNA')
    lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript')
    if raw_type in valid_gene_type: return raw_type;
    if re.match('.*pseudogene',raw_type): return 'pseudogene';
    if raw_type == 'protein_coding': return 'mRNA';
    if raw_type in lncRNA_class: return 'lncRNA';
    return 'other'

def getGeneTypeColor(gType):
    colors = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), (1.0, 0.4980392156862745, 0.054901960784313725), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), (0.8392156862745098, 0.15294117647058825, 0.1568627450980392), (0.5803921568627451, 0.403921568627451, 0.7411764705882353), (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), (0.8901960784313725, 0.4666666666666667, 0.7607843137254902), (0.4980392156862745, 0.4980392156862745, 0.4980392156862745), (0.7372549019607844, 0.7411764705882353, 0.13333333333333333), (0.09019607843137255, 0.7450980392156863, 0.8117647058823529)]
    gTypes = ['mRNA', 'lncRNA', 'pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA']
    if gType == 'other':
        return 'gray'
    i = gTypes.index(gType)
    return colors[i]

def count_valid_ratio(SHAPE):
    ratio_list = []
    for tid in SHAPE:
        total_len = len(SHAPE[tid])
        valid_len = total_len - SHAPE[tid].count('NULL')
        ratio_list.append( 1.0*valid_len/total_len )
    ratio_list.sort()
    return ratio_list

def pie(num_list, textColors=None, textColorFollows=False, explodes=None, colors=None, labels=None, format=None, labeldistance=1.05, fontweight=12):
    import matplotlib.pyplot as plt
    obj = plt.pie(num_list, explode=explodes, colors=colors, labels=labels, shadow=False, autopct=format, labeldistance=labeldistance)
    if textColorFollows:
        for b, t in zip(obj[0], obj[1]):
            t.set_color( b.get_facecolor() )
            t.set_fontsize(12)
    if textColors:
        i = 0
        for b, t in zip(obj[0], obj[1]):
            t.set_color( textColors[i] )
            i += 1
    
    return obj

def rename(rawDict):
    newDict = {}
    for k in rawDict:
        GType = gene_type(k)
        newDict[GType] = newDict.get(GType, 0) + rawDict[k]
    return newDict

def prepare_gene_pie_elements(inDict):
    newDict = rename(inDict)
    
    geneCount = []
    for k in newDict:
        geneCount.append( (k, newDict[k]) )
    geneCount.sort(key=lambda x: x[1], reverse=True)
    gTypes = [ it[0] for it in geneCount ]
    counts = [ it[1] for it in geneCount ]
    ratios = [ 1.0*counts[k]/sum(counts) for k in range(len(counts)) ]
    colors = [ getGeneTypeColor(it) for it in gTypes ]
    
    labels = gTypes[:]
    for i in range(len(ratios)):
        if ratios[i] < 0.01:
            labels[i] = ""
    
    df = []
    for gtype,label,count,color in zip(gTypes, labels, counts, colors):
        df.append( (gtype,label,count,color) )
    df = pd.DataFrame(df, columns=['gtype','label','count','color'])
    df['ratio'] = df['count']/df['count'].sum()
    df['ratioText'] = [ "%.2f%%" % (100*it, ) for it in list(df['ratio']) ]
    
    return df

def classify_trans(SHAPE, Parser):
    TransDict = {}
    BaseDict = {}
    for tid in SHAPE:
        try:
            ft = Parser.getTransFeature(tid)
            gt = ft['gene_type']
            TransDict[ gt ] = TransDict.get(gt, 0) + 1
            valid_num = len(SHAPE[tid]) - SHAPE[tid].count("NULL")
            BaseDict[ gt ] = BaseDict.get(gt, 0) + valid_num
        except KeyError:
            continue
    return TransDict, BaseDict

def cdf(data_list, color='red', topdown=False, label=None, plotMedian=True):
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    
    data_list = np.sort(data_list)
    if topdown:
        p = 1 - 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    else:
        p = 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    plt.plot(data_list, p, color=color, label=label)
    if plotMedian:
        median_x = data_list[ len(data_list)/2 ]
        median_y = p[ len(p)/2 ]
        plt.plot([median_x], [median_y], 'bo')
        plt.axvline(x=median_x, ymin=0, ymax=1, linewidth=2, color='r')

def PlotTransSHAPEStatistics(SHAPE, Parser, outPDF):
    import matplotlib.pyplot as plt
    
    TransDict, BaseDict = classify_trans(SHAPE, Parser)
    df_trans_element = prepare_gene_pie_elements(TransDict)
    df_base_element = prepare_gene_pie_elements(BaseDict)
    
    ratio_list = count_valid_ratio(SHAPE)
    
    print "Plot figures..."
    fig = plt.figure(figsize=(10, 10))
    grids = GridSpec(3, 2)
    
    #### Block 1
    plt.subplot(grids[0, 0], aspect=1)
    objs1 = pie(list(df_trans_element['count']), colors=list(df_trans_element['color']), labels=list(df_trans_element['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Transcripts statistics", fontdict={'fontweight': 'bold'})
    
    #### Block 2
    plt.subplot(grids[0, 1], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_trans_element['color']) ]
    table1 = plt.table(cellText=df_trans_element.loc[:, ('gtype','ratioText')].values, loc='center', cellLoc='left', colLabels=['geneType', 'ratio'], cellColours=cellColors)
    cell_dict = table1.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 3
    plt.subplot(grids[1, 0], aspect=1)
    objs2 = pie(list(df_base_element['count']), colors=list(df_base_element['color']), labels=list(df_base_element['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Bases statistics", fontdict={'fontweight': 'bold'})
    
    #### Block 4
    plt.subplot(grids[1, 1], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_base_element['color']) ]
    table2 = plt.table(cellText=df_base_element.loc[:, ('gtype','ratioText')].values, loc='center', cellLoc='left', colLabels=['geneType', 'ratio'], cellColours=cellColors)
    cell_dict = table2.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 5
    plt.subplot(grids[2, :], aspect=1)
    cdf(ratio_list, color='black', topdown=True, label=None, plotMedian=True)
    plt.xlabel("Cover ratio")
    plt.ylabel("Sorted transcript")
    plt.title("Covered ratio of sorted transcripts")
    
    plt.savefig(outPDF)
    plt.close()











