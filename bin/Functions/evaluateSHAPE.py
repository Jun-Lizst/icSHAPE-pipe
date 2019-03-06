#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import random
import version
import sklearn, sklearn.metrics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Usage = """
evaluateSHAPE - Calculate SHAPE AUC and plot ROC corve with known structure
===========================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i trans_shape.out -s structure.dot -o report.pdf
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input transcript-based SHAPE
  -s                    <String>
                                Input known structure with dot-bracked format
  -o                    <String>
                                Output a PDF report (default: report.pdf)

 More options:
  --step                <Float>
                                0-1, step size (default: 0.01)

\x1b[1mWARNING:\x1b[0m
    typical dot-bracked format:
        >seq1
        ATCGAGTAGCATCGTACGAT
        ....((((....))))....
        >seq2
        GCTGAGTCAGCTAGCTAGCTAAGA
        .((((....).)).((...))...

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inSHAPE': None, 'inDotbracket': None, 'outPDF': 'report.pdf', 'step': 0.01 }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:s:o:', ['step='])
    for op, value in opts:
        if op == '-h':
            print >>sys.stdout, Usage;
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inSHAPE'] = os.path.abspath(value)
        elif op == '-s':
            params['inDotbracket'] = os.path.abspath(value)
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '--step':
            params['step'] = float(value)
            assert 0 < params['step'] < 1
        else:
            print >>sys.stderr, "Error: unrecognized parameter: "+op
            print >>sys.stdout, Usage;
            sys.exit(-1)
    
    # check
    if (not params['inSHAPE']):
          print >>sys.stderr, "Error: Please specify -i"
          print >>sys.stdout, Usage
          sys.exit(-1)
    if (not params['inDotbracket']):
          print >>sys.stderr, "Error: Please specify -s"
          print >>sys.stdout, Usage
          sys.exit(-1)
    
    return params

def Shape_positive_rate(ss_code, shape_list, cutoff):
    Pos_Num = 0
    True_Pos = 0
    False_Pos = 0
    Neg_Num = 0
    for idx, code in enumerate(list(ss_code)):
        if shape_list[idx] != 'NULL':
            if code != ".":
                Pos_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    True_Pos += 1
                else:
                    pass
            else:
                Neg_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    False_Pos += 1
                else:
                    pass
    return 1.0*True_Pos/Pos_Num, 1.0*False_Pos/Neg_Num

def calc_shape_ROC(ss_code, shape_list, step=0.01):
    assert(len(ss_code)==len(shape_list))
    ROC = []
    cutoff = -step
    while cutoff < 1.0 + step:
        TPR, FPR = Shape_positive_rate(ss_code, shape_list, cutoff)
        ROC.append( (FPR, TPR) )
        cutoff += step
    return ROC

def calc_AUC(ROC):
    x = [it[0] for it in ROC]
    y = [it[1] for it in ROC]
    return sklearn.metrics.auc(x, y)

def readDot(inFile):
    structure = {}
    IN = open(inFile)
    line = IN.readline()
    while line:
        if line[0] == '>':
            tid = line[1:].strip().split()[0]
            seq = IN.readline().strip()
            ss = IN.readline().strip()
            structure[tid] = ss
        line = IN.readline()
    IN.close()
    return structure

def readTransSHAPE(file_name):
    transSHAPE = {}
    for line in open(file_name):
        arr = line.strip().split()
        transSHAPE[ arr[0] ] = arr[3:]
    return transSHAPE

def main():
    params = init()
    
    dotbracket = readDot(params['inDotbracket'])
    transSHAPE = readTransSHAPE(params['inSHAPE'])
    
    common_tid = list(set(dotbracket) & set(transSHAPE))
    print "Common transcript in structure file and SHAPE file: ", common_tid
    
    ss_combine = ""
    shape_combine = []
    for tid in common_tid:
        ss = dotbracket[tid]
        shape = transSHAPE[tid]
        if len(shape) != len(ss):
            print >>sys.stderr, "Error: structure length %s != shape length %s in %s" % (len(ss), len(shape), tid)
            exit(-1)
        ss_combine += ss
        shape_combine += shape
    
    ROC = calc_shape_ROC(ss_combine, shape_combine, step=params['step'])
    AUC = calc_AUC(ROC)
    
    fig = plt.figure(1, figsize=(7, 6))
    ax = fig.add_subplot(111)
    
    x = [ i[0] for i in ROC ]
    y = [ i[1] for i in ROC ]
    ax.plot(x, y, '-', color='black')
    ax.plot([0,1], [0,1], '-', color='gray')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC of " + ",".join(common_tid)+" AUC="+str(round(AUC, 2)))
    fig.tight_layout()
    fig.savefig(params['outPDF'])

if __name__ == '__main__':
    main()

