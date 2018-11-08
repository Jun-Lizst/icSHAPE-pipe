#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import commands
import re
import numpy
import version

Usage = """
starbuild - Build a STAR index
=========================================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [--gtf GTF_File -p 1] -i genome.fa -o out_dir
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a genome fasta file
  -o                    <String>
                            Output a path to save index

  More options:
  --gtf                 <String>
                            Provide a GTF file to build index
  -p                    <Int>
                            How many threads to use (default: 1)

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    import getopt
    
    Params = { 'inFile': None, 'outDir': None, 'gtfFile':None, 'threads': 1 }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:p:', ['gtf='])
    
    for op, value in opts:
        if op == '-h':
            print Usage
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            Params['outDir'] = os.path.abspath(value.rstrip('/'))
        
        elif op == '--gtf':
            Params['gtfFile'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        
        else:
            print >>sys.stderr, "parameter Error: unrecognized parameter: "+op
            print Usage
            sys.exit(-1)
    
    if Params['inFile'] == None:
        print >>sys.stderr, "Error: please specify -i"
        print Usage
        exit(-1)
    if Params['outDir'] == None:
        print >>sys.stderr, "Error: please specify -o"
        print Usage
        exit(-1)
    
    return Params

def count_fasta(inFasta):
    ref_num = 0
    base_num = 0
    for line in open(inFasta):
        if line[0] == '>':
            ref_num += 1
        else:
            base_num += len(line) - 1
    return ref_num, base_num


def main():
    params = init()
    CMD = "STAR --runMode genomeGenerate --genomeFastaFiles %s --genomeDir %s --runThreadN %s" % (params['inFile'], params['outDir'], params['threads'])
    if params['gtfFile']:
        CMD += " --sjdbGTFfile " + params['gtfFile']
    
    ref_num, base_num = count_fasta(params['inFile'])
    genomeSAindexNbases = min( 14, numpy.log2(base_num)/2-1 )
    if genomeSAindexNbases < 14:
        CMD += " --genomeSAindexNbases %.1f" % (genomeSAindexNbases, )
    
    if ref_num > 5000:
        genomeChrBinNbits = min( 18, numpy.log2(1.0*base_num/ref_num) )
        if genomeChrBinNbits < 18:
            CMD += " --genomeChrBinNbits %.1f" % (genomeChrBinNbits, )
    
    print "Start to build STAR index:\n\t%s" % (CMD, )
    os.system(CMD)

if __name__ == "__main__":
    main()
