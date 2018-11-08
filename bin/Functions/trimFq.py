#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import commands
import version

Usage = """
trimFq - Trim 5' barcode and 3' adaptor
=========================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1 -m 25] -i inFastq -o outFastq -l trim_leading_len -a adaptor
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a fastq file
  -o                    <String>
                            Output a processed fastq file
  -l                    <Int>
                            How many base to trim 5' of reads (typical icSHAPE: 13)
  -a                    <String>
                            Input a adaptor fasta file

  More options:
  -p                    <Int>
                            How many threads to use (default: 1)
  -m                    <Int>
                            Minimun length to preserve (default: 25)

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

dirname = os.path.dirname(os.path.abspath(__file__))
trimmomatic = os.path.join(dirname, 'trimmomatic-0.38.jar')

def init():
    import getopt
    
    Params = { 'inFastq': None, 'outFastq': None, 'leading': 13, 'adator': None, 'threads': 1, 'minLen': 25 }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:l:a:p:m:')
    
    for op, value in opts:
        if op == '-h':
            print Usage
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFastq'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFastq'] = os.path.abspath(value)
        elif op == '-l':
            Params['leading'] = int(value)
        elif op == '-a':
            Params['adator'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        elif op == '-m':
            Params['minLen'] = int(value)
        
        else:
            print >>sys.stderr, "parameter Error: unrecognized parameter: "+op
            print Usage
            sys.exit(-1)
    
    if not Params['inFastq']:
        print >>sys.stderr, "Error: please specify -i"
        print Usage
        exit(-1)
    if not Params['outFastq']:
        print >>sys.stderr, "Error: please specify -o"
        print Usage
        exit(-1)
    if not Params['adator']:
        print >>sys.stderr, "Error: please specify -a"
        print Usage
        exit(-1)
    
    return Params


def main():
    params = init()
    CMD = "java -mx256m -jar %s SE -threads %s -phred33 %s %s HEADCROP:%s ILLUMINACLIP:%s:2:30:4 TRAILING:20 MINLEN:%s" % (trimmomatic, params['threads'], params['inFastq'], params['outFastq'], params['leading'], params['adator'], params['minLen'])
    print "Start to build STAR index:\n\t%s" % (CMD, )
    os.system(CMD)

if __name__ == "__main__":
    main()





