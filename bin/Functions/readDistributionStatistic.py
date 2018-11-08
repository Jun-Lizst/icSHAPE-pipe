#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import version
import matplotlib.pyplot as plt

Usage = """
readDistributionStatistic - Count fastq or sam to statistic reads distribution
==============================================================================
\x1b[1mUSAGE:\x1b[0m
  %s [--labels D1,D2,N1,N2] -1 raw_fq.fastq -2 read_collapse.fastq -3 trimmed.fastq -4 rRNA.fastq -5 mapGenome.sam -o report.pdf

\x1b[1mHELP:\x1b[0m
  -1                <String,String...>
                            Raw unprocessed fastq file
  -2                <String,String...>
                            Collapsed fastq file
  -3                <String,String...>
                            Trimmed fastq file
  -4                <String,String...>
                            rRNA removed fastq file
  -5                <String,String...>
                            Sam file after mapping to genome
  -o                <String>
                            Output a PDF report file (default: report.pdf)

 More options:
  --labels          <String,String...>
                            Label for each sample

\x1b[1mWARNING:\x1b[0m
    1. Multiple samples can be concated by comma;
    2. Compressed files (.gz and bam) are permitted.

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan


""" % (sys.argv[0], version.Version)


def init():
    params = { '1': [], '2': [], '3': [], '4': [], '5': [], 'outPDF': 'report.pdf', 'labels': [] }
    opts, args = getopt.getopt(sys.argv[1:], 'h1:2:3:4:5:o:', ['labels='])
    for op, value in opts:
        if op == '-h':
            print >>sys.stdout, Usage;
            sys.exit(-1)
        # Basic Parameters
        elif op == '-1':
            for file in value.strip().split(','):
                params['1'].append( os.path.abspath(file) )
        elif op == '-2':
            for file in value.strip().split(','):
                params['2'].append( os.path.abspath(file) )
        elif op == '-3':
            for file in value.strip().split(','):
                params['3'].append( os.path.abspath(file) )
        elif op == '-4':
            for file in value.strip().split(','):
                params['4'].append( os.path.abspath(file) )
        elif op == '-5':
            for file in value.strip().split(','):
                params['5'].append( os.path.abspath(file) )
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '--labels':
            for label in value.strip().split(','):
                params['labels'].append(label)
        else:
            print >>sys.stderr, "Error: unrecognized parameter: "+op
            print >>sys.stdout, Usage;
            sys.exit(-1)
    
    # check
    sn = len(params['1'])
    if (sn==0) or (sn!=len(params['2'])) or (sn!=len(params['3'])) or (sn!=len(params['4'])) or (sn!=len(params['5'])):
          print >>sys.stderr, "Error: file nums or not equal"
          print >>sys.stdout, Usage
          sys.exit(-1)
    if len(params['labels']) == 0:
        params['labels'] = [ str(i) for i in range(1, sn+1) ]
    elif sn != len(params['labels']):
        print >>sys.stderr, "Error: --labels different length"
        print >>sys.stdout, Usage
        sys.exit(-1)
    
    return params

def count_fq(inFile):
    import commands
    
    if inFile.endswith('.gz'):
        count = commands.getoutput("gzip -d %s -c | wc -l" % (inFile, ))
        count = int(count)
    else:
        count = commands.getoutput("wc -l %s | cut -f 1 -d \" \"" % (inFile, ))
        count = int(count)
    
    assert count%4 == 0
    return count/4

def count_sam(inFile):
    import commands
    
    CMD = "samtools view %s | awk 'BEGIN{map=0;unmap=0;}{if($3==\"*\"){unmap+=1}else{map+=1}}END{print map\"\t\"unmap}'"
    
    mapcount, unmapcount = commands.getoutput(CMD % (inFile, )).strip().split('\t')
    return int(mapcount), int(unmapcount)


def stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors=None ):
    """
    stackedBars = [ [29, 10, 21], [24, 11, 33] ]
    stackedLabels = ['stack1', 'stack2', 'stack3'] 
    barLabels = ['bar1', 'bar2']
    stackedColors = ['red', 'blue', 'green']
    stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors)
    """
    
    import matplotlib.pyplot as plt
    
    def getStackedBarList(stackedBars, stackedLabels, barLabels):
        barList = []
        for bar,bar_label in zip(stackedBars,barLabels):
            for bar_item,stack_label in zip(bar,stackedLabels):
                barList.append( (bar_item, stack_label, bar_label) )
        return barList
    
    assert len(stackedBars) == len(barLabels)
    assert  len(stackedBars[0]) == len(stackedLabels)
    
    if stackedColors:
        assert len(stackedLabels) == len(stackedColors)
    else:
        stackedColors = sns.color_palette("hls", len(stackedLabels))
    
    stackedBarList = getStackedBarList(stackedBars, stackedLabels, barLabels)
    
    last_y = [0]*len(barLabels)
    for i,stack_label in enumerate( stackedLabels ):
        sub_dict = {it[2]:it[0] for it in stackedBarList if it[1]==stack_label}
        
        y = []
        for bar_label in barLabels:
            y.append(sub_dict[bar_label])
        
        plt.bar(range(1,len(barLabels)+1), y, color=stackedColors[i], bottom=last_y, label=stack_label)
        last_y = [ y_i+y_j for y_i,y_j in zip(y, last_y) ]
    
    plt.legend(frameon=False)
    plt.xticks( range(1,len(barLabels)+1), barLabels )
    plt.xticks(rotation=90)


def main():
    params = init()
    sn = len(params['1'])
    stackedBars = []
    for i in range(sn):
        print "Start to statistic sample "+str(i)+"..."
        
        raw_count = count_fq(params['1'][i])
        uniq_count = count_fq(params['2'][i])
        trim_count = count_fq(params['3'][i])
        clean_count = count_fq(params['4'][i])
        map_count,unmap_count = count_sam(params['5'][i])
        
        dup_num = raw_count - uniq_count
        trim_num = uniq_count - trim_count
        rRNA_num = trim_count - clean_count
        
        if raw_count != dup_num+trim_num+rRNA_num+map_count+unmap_count:
            print >>sys.stderr, "Error: raw_count != dup_num+trim_num+rRNA_num+map_count+unmap_count"
            exit(-1)
        
        stackedBars.append( (dup_num, trim_num, rRNA_num, map_count, unmap_count) )
    
    print "Start to plot..."
    stackedLabels = ['Duplicates', 'Trimmed', 'rRNA_tRNA_mtRNA', 'Map_genome', 'Unmap'] 
    barLabels = params['labels']
    stackedColors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974']
    plt.figure(figsize=(8, 10))
    stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors)
    plt.xlabel("Sample")
    plt.ylabel("Number of reads")
    plt.tight_layout()
    plt.savefig(params['outPDF'])


if __name__ == '__main__':
    main()


