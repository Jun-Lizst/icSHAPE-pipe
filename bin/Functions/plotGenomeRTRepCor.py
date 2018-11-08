#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import random
import version


Usage = """
plotGenomeRTRepCor - Calculate replicate correlation for genome RT
==================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i inputFile -o report.pdf --col1 4 --col2 6
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                input countRT file record the RT (prduced by countRT)
  -o                    <String>
                                output a PDF report (default: report.pdf)
  --col1                <Int>
                                RT column number of replicate 1 (default: 4)
  --col2                <Int>
                                RT column number of replicate 2 (default: 6)

  More options:
  --minBD               <Int>
                                minimun basedensity (default: 100)
  --winSize             <Int>
                                window size for each replicate calculation (default: 100)

\x1b[1mWARNING:\x1b[0m
    Basedensity is appended after the RT column

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inFile': None, 'outPDF': 'report.pdf', 'col1': 4, 'col2': 6, 'minBD': 100, 'winSize': 100 }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['col1=', 'col2=', 'minBD=', 'winSize='])
    for op, value in opts:
        if op == '-h':
            print >>sys.stdout, Usage;
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '--col1':
            params['col1'] = int(value)
        elif op == '--col2':
            params['col2'] = int(value)
        elif op == '--minBD':
            params['minBD'] = int(value)
        elif op == '--winSize':
            params['winSize'] = int(value)
        else:
            print >>sys.stderr, "Error: unrecognized parameter: "+op
            print >>sys.stdout, Usage;
            sys.exit(-1)
    
    # check
    if (not params['inFile']):
          print >>sys.stderr, "Error: Please specify -i"
          print >>sys.stdout, Usage
          sys.exit(-1)
    
    return params


def calcRTReplicateCorrelation(inFile, col1, col2, minBD=100, windowsize=100):
    import scipy
    import scipy.stats
    
    rt1 = []
    rt2 = []
    true_cor = []
    flip1_cor = []
    rand_cor = []
    
    for line in open(inFile):
        data = line.strip().split()
        r1, r2 = int(data[col1-1]), int(data[col2-1])
        b1, b2 = int(data[col1]), int(data[col2])
        if b1 < minBD or b2 < minBD:
            continue
        
        rt1.append(r1)
        rt2.append(r2)
        if len(rt1) == windowsize:
            p, v = scipy.stats.pearsonr(rt1, rt2)
            if v < 0.05:
                true_cor.append( p )
            
            p, v = scipy.stats.pearsonr(rt1, [rt2[-1]]+rt2[:-1])
            if v < 0.05:
                flip1_cor.append( p )
            
            shuffle_rt2 = rt2[:]
            random.shuffle(shuffle_rt2)
            p, v = scipy.stats.pearsonr(rt1, shuffle_rt2)
            if v < 0.05:
                rand_cor.append( p )
            
            rt1 = []
            rt2 = []
    
    return true_cor, flip1_cor, rand_cor


def boxplot(data_list, ax, width=0.4, labels=None, title=None):
    """
    Example:
        fig = plt.figure(1, figsize=(5, 6))
        ax = fig.add_subplot(111)
        boxplot(data_list, ax)
        fig.show()
    """
    obj = ax.boxplot(x=data_list, showfliers=False, patch_artist=True, widths = width)
    for box in obj['boxes']:
        box.set( color='#7570b3', linewidth=0.5)
        box.set( facecolor = '#1b9e77' )
    
    for whisker in obj['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    
    for cap in obj['caps']:
        cap.set(color='#7570b3', linewidth=2)
    
    for median in obj['medians']:
        median.set(color='#b2df8a', linewidth=2)
    
    for flier in obj['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    
    if labels:
        ax.set_xticklabels(labels)
    
    if title:
        ax.set_title(title).set_weight("bold")
    
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    return obj


def main():
    params = init()
    import matplotlib.pyplot as plt
    true_cor, flip1_cor, rand_cor = calcRTReplicateCorrelation(inFile=params['inFile'], col1=params['col1'], col2=params['col2'], minBD=params['minBD'], windowsize=params['winSize'])
    
    fig = plt.figure(1, figsize=(5, 6))
    ax = fig.add_subplot(111)
    boxplot([true_cor, flip1_cor, rand_cor], ax, width=0.4, labels=['Pearson\ncorrelation', 'Flip-1bp\ncorrelation', 'Shuffled\ncorrelation'])
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("RT pearson correlation for replicates").set_weight("bold")
    ax.set_ylabel("Pearson correlation for each window").set_weight("bold")
    fig.savefig(params['outPDF'])

if __name__ == '__main__':
    main()





"""

gTab1 = "/150T/zhangqf/lipan/icSHAPE_Pipeline/new_pipeline/7.replicateControl/shape_1.gTab"
gTab2 = "/150T/zhangqf/lipan/icSHAPE_Pipeline/new_pipeline/7.replicateControl/shape_2.gTab"
outFile = "/150T/zhangqf/lipan/icSHAPE_Pipeline/new_pipeline/7.replicateControl/combine.txt"

combine_gTab(gTab1, gTab2, outFile)



true_cor_50, rand_cor_50 = calcSHAPEReplicateCorrelation(outFile, windowsize=50); print len(true_cor_50),len(rand_cor_50)
true_cor_100, rand_cor_100 = calcSHAPEReplicateCorrelation(outFile, windowsize=100); print len(true_cor_100),len(rand_cor_100)
true_cor_200, rand_cor_200 = calcSHAPEReplicateCorrelation(outFile, windowsize=200); print len(true_cor_200),len(rand_cor_200)

print numpy.median(true_cor_50), numpy.median(rand_cor_50)
print numpy.median(true_cor_100), numpy.median(rand_cor_100)
print numpy.median(true_cor_200), numpy.median(rand_cor_200)



inFile = "/150T/zhangqf/lipan/icSHAPE_Pipeline/new_pipeline/7.replicateControl/countRT.txt"
true_cor_50, rand_cor_50 = calcRTReplicateCorrelation(inFile, col1=4, col2=6, windowsize=50); print len(true_cor_50),len(rand_cor_50)
true_cor_100, rand_cor_100 = calcRTReplicateCorrelation(inFile, col1=4, col2=6, windowsize=100); print len(true_cor_100),len(rand_cor_100)
true_cor_200, rand_cor_200 = calcRTReplicateCorrelation(inFile, col1=4, col2=6, windowsize=200); print len(true_cor_200),len(rand_cor_200)


print numpy.median(true_cor_50), numpy.median(rand_cor_50)
print numpy.median(true_cor_100), numpy.median(rand_cor_100)
print numpy.median(true_cor_200), numpy.median(rand_cor_200)

"""