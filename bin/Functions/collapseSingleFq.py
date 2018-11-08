#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import commands
import gzip
import version

Usage = """
collaseSingleFq - Collapse a fastq file
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i input_fastq -o output_fastq --mode mode --fasta output_fasta
\x1b[1mHELP:\x1b[0m
    --mode:             <append/new>. 
                            append will lump reads into files in the outputdirectory with names like LIB_DMSO1.fastq, LIB_DMSO2.fastq...
                            new will create separate files for each input file, i.e., names like LIB_DMSO1_dataset1.fastq, LIB_DMSO2_dataset1.fastq...
    --fasta:            <String>
                            output corresponding fasta file

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    import getopt
    
    Params = { 'mode': 'new', 'outFasta': None }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['mode=', 'fasta='])
    
    for op, value in opts:
        if op == '-h':
            print Usage
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFile'] = os.path.abspath(value)
        elif op == '--mode':
            Params['mode'] = value
        elif op == '--fasta':
            Params['outFasta'] = os.path.abspath(value)
        else:
            print >>sys.stderr, "parameter Error: unrecognized parameter: "+op
            print Usage
            sys.exit(-1)
    
    if 'inFile' not in Params:
        print >>sys.stderr, "Error: specify -i"
        print Usage
        exit(-1)
    if 'outFile' not in Params:
        print >>sys.stderr, "Error: specify -o"
        print Usage
        exit(-1)
    
    return Params

def main():
    Params = init()
    
    print >>sys.stderr, "Collapsing file %s...\n\t" % (Params['inFile'], )
    
    totalCount, representativeCount, collapseRatio, fastaFile = collapse( Params['inFile'], Params['outFile'], Params['mode'], Params['outFasta'] )
    
    print >>sys.stderr, "Read collapse successful! Total count: %s, unique count: %s, unique ratio: %s. \nRefer to the collapsed fasta file: %s" % (totalCount, representativeCount, collapseRatio, fastaFile)
    print >>sys.stderr, "Collapsing file %s finished." % (Params['inFile'], ) 
    
    return 0

def collapse(inFile, outFile, mode="new", fastaFile=None, compress=False):
    
    fileName, fileDir, fileSuffix = fileparse(outFile)
    fastaFile = fastaFile if fastaFile else fileDir+"/"+fileName+".fa"
    
    CL = gzip.open(inFile, 'r') if inFile.endswith('.gz') else open(inFile, 'r')
    
    if mode == "new":
        OUT1 = gzip.open(outFile, 'w') if outFile.endswith('.gz') else open(outFile, 'w')
        FA = gzip.open(fastaFile, 'w') if fastaFile.endswith('.gz') else open(fastaFile, 'w')
    elif mode == 'append':
        OUT1 = gzip.open(outFile, 'a') if outFile.endswith('.gz') else open(outFile, 'a')
        FA = gzip.open(fastaFile, 'a') if fastaFile.endswith('.gz') else open(fastaFile, 'a')
    
    readCount = {}
    totalCount = 0
    representativeCount = 0
    key = CL.readline()
    while key:
        seq = CL.readline()
        tag = CL.readline()
        quality = CL.readline()
        
        totalCount += 1
        if totalCount % 1000000 == 0:
            print >>sys.stderr, totalCount
        
        try:
            readCount[seq] += 1
        except KeyError:
            readCount[seq] = 1
            OUT1.writelines( key+seq+tag+quality )
            representativeCount += 1
        
        key = CL.readline()
    CL.close()
    OUT1.close()

    FA.writelines("## -------------------------------------------------------------------------\n");
    FA.writelines("## fasta sequence frequencies from file %s.\n" % (inFile, ))
    FA.writelines("## ------------------\n")
    sorted_key = sorted( readCount.keys(), key=lambda x: readCount[x], reverse=True )
    for seq_key in sorted_key:
        FA.writelines( "%s\t%s\n" % (seq_key.strip(), readCount[seq_key]) )
    FA.close()
    
    collapseRatio = "%.2f" % (1.0*representativeCount/totalCount, )
    return totalCount, representativeCount, collapseRatio, fastaFile


def fileparse(fillFilePath):
    fileDir = '/'.join(fillFilePath.split('/')[:-1])
    fileName = fillFilePath.split('/')[-1]
    pureFileName = '.'.join(fileName.split('.')[:-1])
    fileSuffix = fileName.split('.')[-1]
    return pureFileName, fileDir, fileSuffix

if __name__ == "__main__":
    main()

