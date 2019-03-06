#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os, commands, sys, pysam, numpy, getopt, random, re
import GAP
from matplotlib.gridspec import GridSpec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import version

from StatisticPlot import *

Usage = """
samStatistic - Report the sam or bam
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-s 0.1] -i input.sam -g *.genomeCoor.bed -o report.pdf -t report.txt
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input a sam or bam file, each read must be output for one time
  -o                    <String>
                                Output statistics to PDF file (default: report.pdf)
  -t                    <String>
                                Output statistics to text file (default: report.txt)
  --fast                <None>
                                Fast mode, sample reads to analysis (default: all reads will be analized)
  -g                    <String>
                                A genome-coor based annotation file: hg38.genomeCoor.bed (generated by parseGTF)

\x1b[1mWARNING\x1b[0m
    1. Each mapped read should be output for one time (--outSAMmultNmax 1);
    2. Each mapped read should contain MD tag (such as MD:Z:35);
    3. Unmapped reads should be append in sam/bam file (--outSAMunmapped Within);
    4. Annotation come from GAP.

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    params = { 'inFile': None, 'outPDF': 'report.pdf', 'outTXT': 'report.txt', 'fast': False, 'annotation': None }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:g:t:', ['fast'])
    for op, value in opts:
        if op == '-h':
            print >>sys.stdout, Usage;
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '-t':
            params['outTXT'] = os.path.abspath(value)
        elif op == '--fast':
            params['fast'] = True
        elif op == '-g':
            params['annotation'] = os.path.abspath(value)
        
        else:
            print >>sys.stderr, "Error: unrecognized parameter: "+op
            print >>sys.stdout, Usage;
            sys.exit(-1)
    # check
    if not params['inFile'] or not params['annotation']:
        print >>sys.stderr, "Error: Please specify -i -g"
        print >>sys.stdout, Usage
        sys.exit(-1)
    
    return params

def findOverlapGeneTypes(Pos, tIDList, Parser):
    geneTypes = set()
    for tid in tIDList:
        tFeature = Parser.getTransFeature(tid)
        tStart = tFeature['start']
        tEnd = tFeature['end']
        if tStart <= Pos <= tEnd:
            geneTypes.add( tFeature['gene_type'] )
    return geneTypes


def estimate_sample_ratio(inBam):
    CMD = "samtools view %s | wc -l" % (inBam, )
    count = int(commands.getoutput(CMD))
    if count < 1000000:
        return None
    else:
        return 1.0 * 1000000 / count

def get_mute_list(MDCode):
    mute_list = []
    
    match_list = re.split("[ATGC^N]", MDCode)[:-1]
    ext_len = 0
    for cur_len in match_list:
        cur_len = 1 if cur_len == '' else int(cur_len)
        ext_len += cur_len
        mute_list.append(ext_len)
        ext_len += 1
    
    return mute_list

def evaluate_sam(inSam, Parser, sample=None):
    
    geneParser = Parser.getGeneParser()
    
    ## Mapped reads classification
    TotalCount = 0
    Unmap = 0
    Exon = { 'mRNA': 0 }
    Intron = { 'mRNA': 0 }
    Intergenic = 0
    Unannot = 0
    
    ## Mutation profile
    mutTotalCount = 0
    mutProf_leading = [0] * 15
    mutProf_tailing = [0] * 15
    
    ## Open file
    if inSam.endswith('sam'):
        IN = pysam.AlignmentFile(inSam, "r")
    elif inSam.endswith('bam'):
        IN = pysam.AlignmentFile(inSam, "rb")
    else:
        print >>sys.stderr, "Error: input must be a sam or bam file"
        exit(-1)
    
    last_query_id = ""
    while 1:
        try: alignSeg = IN.next()
        except StopIteration: break
        
        cur_query_name = alignSeg.query_name
        if cur_query_name == last_query_id:
            print >>sys.stderr, "Error: Each reads should be output for one time (--outSAMmultNmax 1)"
            exit(-1)
        last_query_id = cur_query_name
        
        TotalCount += 1
        if TotalCount % 100000 == 0:
            print "lines ", TotalCount
        
        if sample:
            if random.random() > sample:
                continue
        
        if alignSeg.is_unmapped:
            Unmap += 1
            continue
        
        mutTotalCount += 1
        MDCode = dict(alignSeg.tags)['MD']
        muteList = get_mute_list(MDCode)
        seq_len = len(alignSeg.seq)
        for coor in muteList:
            if coor < len(mutProf_leading):
                mutProf_leading[coor] += 1
            else:
                index = (coor+1)-(seq_len-15)
                if 0 <= index < len(mutProf_tailing):
                    mutProf_tailing[index] += 1
        
        chrID = alignSeg.reference_name
        chrPos = alignSeg.pos + 1
        strand = '-' if alignSeg.is_reverse else '+'
        
        try:
            geneList = Parser.genomeCoor2geneCoor(chrID, chrPos, chrPos, strand)
        except KeyError:
            Unannot += 1
            continue
        
        if len(geneList)==0:
            Intergenic += 1
        else:
            transList = Parser.genomeCoor2transCoor(chrID, chrPos, chrPos, strand)
            
            gID = geneList[0][3]
            gTransList = geneParser[gID]['transcript']
            
            gTypes = geneParser[gID]['gene_type']
            if ('protein_coding' in gTypes) or ('mRNA' in gTypes):
                gType = 'mRNA'
            else:
                gType = None
            
            if len(transList)==0:
                if gType != 'mRNA':
                    geneTypes = findOverlapGeneTypes(chrPos, gTransList, Parser)
                    if len(geneTypes) == 0: Unannot += 1; continue
                    add_ratio = 1.0/len(geneTypes)
                    for gt in geneTypes:
                        Intron[ gt ] = Intron.get(gt, 0) + add_ratio
                else:
                    Intron[ 'mRNA' ] += 1
            else:
                if gType != 'mRNA':
                    geneTypes = findOverlapGeneTypes(chrPos, gTransList, Parser)
                    if len(geneTypes) == 0: Unannot += 1; continue
                    add_ratio = 1.0/len(geneTypes)
                    for gt in geneTypes:
                        Exon[ gt ] = Exon.get(gt, 0) + add_ratio
                else:
                    Exon[ 'mRNA' ] += 1
    
    for i in range(len(mutProf_leading)):
        mutProf_leading[i] = round(1.0*mutProf_leading[i]/mutTotalCount, 3)
    
    for i in range(len(mutProf_tailing)):
        mutProf_tailing[i] = round(1.0*mutProf_tailing[i]/mutTotalCount, 3)
    
    return TotalCount, Unmap, Exon, Intron, Intergenic, Unannot, mutProf_leading, mutProf_tailing


def build_sorted_list(inDict):
    myList = []
    for k in inDict:
        myList.append((k, inDict[k]))
    
    myList.sort(key=lambda x: x[1], reverse=True)
    return myList

def writeTextReport(outFile, TotalCount, Unmap, Exon, Intron, Intergenic, Unannot, mutProf_leading, mutProf_tailing):
    Map_list = build_sorted_list({'unmap':Unmap, 'exon':sum(Exon.values()), 'intron':sum(Intron.values()), 'intergenic':Intergenic})
    df_exon = prepare_gene_pie_elements(Exon)
    df_intron = prepare_gene_pie_elements(Intron)
    
    AnaCount = Unmap + sum(Exon.values()) + sum(Intron.values()) + Intergenic + Unannot
    AnnotCount = AnaCount - Unannot
    
    list_exon = df_exon.loc[:, ('gtype', 'count', 'ratio')].values
    list_intron = df_intron.loc[:, ('gtype', 'count', 'ratio')].values
    
    OUT = open(outFile, 'w')
    print >>OUT, "Total reads: %s" % (TotalCount, )
    print >>OUT, "Sampled analyzed reads: %s, ratio: %.2f%%" % (AnaCount, 100.0*AnaCount/TotalCount)
    
    print >>OUT, "Mapped reads: %s, ratio: %.2f%%" % (AnaCount-Unmap, 100-100.0*Unmap/AnaCount)
    print >>OUT, "\tReads map to exons: %s, ratio: %.2f%%" % (sum(Exon.values()), 100.0*sum(Exon.values())/AnnotCount)
    for line in list_exon:
        print >>OUT, "\t\t%s: %s, ratio: %.2f%%" % (line[0], int(line[1]), 100*line[2])
    
    print >>OUT, "\tReads map to introns: %s, ratio: %.2f%%" % (sum(Intron.values()), 100.0*sum(Intron.values())/AnnotCount)
    for line in list_intron:
        print >>OUT, "\t\t%s: %s, ratio: %.2f%%" % (line[0], int(line[1]), 100*line[2])
    
    print >>OUT, "\tReads map to intergenic regions: %s, ratio: %.2f%%" % (Intergenic, 100.0*Intergenic/AnnotCount)
    
    print >>OUT, "Mutate profile leading 15nt: " + ",".join([str(it) for it in mutProf_leading])
    print >>OUT, "Mutate profile tailing 15nt: " + ",".join([str(it) for it in mutProf_tailing])
    
    OUT.close()


def main():
    param = init()
    
    print "Start to read annotation file..."
    Parser = GAP.init(param['annotation'])
    
    ratio = None
    if param['fast']:
        print "Start to estimate sample ratio..."
        ratio = estimate_sample_ratio(param['inFile'])
        print "Fast mode, sample ratio: %s" % (ratio, )
    
    TotalCount, Unmap, Exon, Intron, Intergenic, Unannot, mutProf_leading, mutProf_tailing = evaluate_sam(param['inFile'], Parser, sample=ratio)
    
    print "Write text report..."
    writeTextReport(param['outTXT'], TotalCount, Unmap, Exon, Intron, Intergenic, Unannot, mutProf_leading, mutProf_tailing)
    
    Map_list = build_sorted_list({'unmap':Unmap, 'exon':sum(Exon.values()), 'intron':sum(Intron.values()), 'intergenic':Intergenic})
    df_exon = prepare_gene_pie_elements(Exon)
    df_intron = prepare_gene_pie_elements(Intron)
    
    ####################
    ###   Plot
    ####################
    
    print "Plot figures..."
    fig = plt.figure(figsize=(10, 15))
    grids = GridSpec(4, 2)
    
    #### Block 1
    plt.subplot(grids[0, 0], aspect=1)
    objs1 = pie(list(df_exon['count']), colors=list(df_exon['color']), labels=list(df_exon['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Mapped to exons", fontdict={'fontweight': 'bold'})
    
    #### Block 2
    plt.subplot(grids[0, 1], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_exon['color']) ]
    table1 = plt.table(cellText=df_exon.loc[:, ('gtype','ratioText')].values, loc='center', cellLoc='left', colLabels=['geneType', 'ratio'], cellColours=cellColors)
    cell_dict = table1.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 3
    plt.subplot(grids[1, 0], aspect=1)
    objs2 = pie(list(df_intron['count']), colors=list(df_intron['color']), labels=list(df_intron['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Mapped to introns", fontdict={'fontweight': 'bold'})
    
    #### Block 4
    plt.subplot(grids[1, 1], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_intron['color']) ]
    table2 = plt.table(cellText=df_intron.loc[:, ('gtype','ratioText')].values, loc='center', cellLoc='left', colLabels=['geneType', 'ratio'], cellColours=cellColors)
    cell_dict = table2.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 5
    plt.subplot(grids[2, 0], aspect=1)
    labels3 = [ it[0] for it in Map_list ]
    values3 = [ it[1] for it in Map_list ]
    faceMap = { 'unmap':'#B4B4B5', 'exon':'#FCB9A9', 'intron':'#9393F9', 'intergenic':'#9BF9B6' }
    textMap = { 'unmap':'#4B4B4C', 'exon':'#D8654E', 'intron':'#4545CE', 'intergenic':'#37BC59' }
    faceColor = [faceMap[k] for k in labels3]
    labelColor = [textMap[k] for k in labels3]
    objs3 = pie(values3, labels=labels3, colors=faceColor, textColorFollows=True, textColors=labelColor, labeldistance=1.03)
    plt.title("Mapped to genome", fontdict={'fontweight': 'bold'})
    
    #### Block 6
    celltext = []
    for gType, ratio in zip(labels3, values3):
        ratioText = "%.2f%%" % ( 100.0*ratio/sum(values3), )
        celltext.append( (gType, ratioText) )
    
    plt.subplot(grids[2, 1], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in faceColor ]
    table3 = plt.table(cellText=celltext, loc='center', cellLoc='left', colLabels=['Type', 'Ratio'], cellColours=cellColors)
    cell_dict = table3.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 7
    ax4 = fig.add_subplot(grids[3, :])
    mutProf = mutProf_leading + [0, 0, 0] + mutProf_tailing
    xticks = [str(i) for i in range(1,16)] + ["", "", ""] + [str(i) for i in range(-15,0)]
    obj = ax4.bar(range(len(mutProf)), mutProf, color='#B4B4B5')
    ax4.title.set_text("Nucleotide mutation profile")
    ax4.title.set_size('large')
    ax4.title.set_weight('bold')
    ax4.set_xticks(range(len(mutProf)))
    ax4.set_xticklabels(xticks)
    ax4.set_ylabel("Mutation ratio")
    ax4.set_xlabel("Nucleotide position in reads")
    
    plt.savefig(param['outPDF'])
    plt.close()

if __name__ == '__main__':
    main()



"""

param = {'inFile': '/150T/zhangqf/lipan/icSHAPE_Pipeline/test/test.Aligned.out.bam', 
         'outPDF': 'figs/result.pdf', 
         'outTXT': 'figs/result.txt', 
         'fast': True,
         'annotation': '/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed'}

"""
