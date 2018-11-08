#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import commands
import version

Usage = """
mapGenome - Map reads to genome with STAR
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1 --maxMMap 1 --maxMisMatch 2 --noMut5] -i inFastq -o outprefix -x index
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a fastq file
  -o                    <String>
                            Prefix (with full path) of all outputfiles (like D1, D2...)
                            Final bam file will be Prefix.sorted.bam
  -x                    <String>
                            Input a index to map

  More options:
  -p                    <Int>
                            How many threads to use (default: 1)
  --maxMMap             <Int>
                            Maximun multiple map to be allowed (default: 1, unique map)
  --maxMisMatch         <Int>
                            Maximun mismatch to be allowed (default: 2)
  --noMut5              <None>
                            Remove reads with mutation at the first base in 5' (such as MD:Z:0A)

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    import getopt
    
    Params = { 'inFastq': None, 'outPrefix': None, 'index': None, 'threads': 1, 'maxMMap': 1, 'maxMisMatch': 2, 'noMut5':False }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:x:p:', ['maxMMap=', 'maxMisMatch=', 'noMut5'])
    
    for op, value in opts:
        if op == '-h':
            print Usage
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFastq'] = os.path.abspath(value)
        elif op == '-o':
            Params['outPrefix'] = os.path.abspath(value)
        elif op == '-x':
            Params['index'] = os.path.abspath(value)
        elif op == '-p':
            Params['adator'] = os.path.abspath(value)
        
        elif op == '--maxMMap':
            Params['maxMMap'] = int(value)
        elif op == '--maxMisMatch':
            Params['maxMisMatch'] = int(value)
        elif op == '--noMut5':
            Params['noMut5'] = True
        
        else:
            print >>sys.stderr, "parameter Error: unrecognized parameter: "+op
            print Usage
            sys.exit(-1)
    
    if not Params['inFastq']:
        print >>sys.stderr, "Error: please specify -i"
        print Usage
        exit(-1)
    if not Params['outPrefix']:
        print >>sys.stderr, "Error: please specify -o"
        print Usage
        exit(-1)
    if not Params['index']:
        print >>sys.stderr, "Error: please specify -x"
        print Usage
        exit(-1)
    
    return Params


def main():
    CMD_1 = "STAR --readFilesIn %s \
        --outFileNamePrefix %s. \
        --genomeDir %s \
        --runThreadN %s \
        --genomeLoad NoSharedMemory \
        --runMode alignReads \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax 1 \
        --outFilterMultimapNmax %s \
        --outFilterMismatchNmax %s \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField intronMotif \
        --outSJfilterOverhangMin 30 12 12 12 \
        --alignEndsType EndToEnd \
        --outSAMattributes All \
        --outSAMunmapped Within \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --outStd BAM_Unsorted > %s"
    
    CMD_sort_1 = "samtools sort %s -o %s"
    CMD_sort_2 = "samtools view -h %s | grep -v \"MD:Z:0\"| samtools view -bh | samtools sort -o %s"
    
    params = init()
    unsorted_bam = params['outPrefix'] + ".unsorted.bam"
    sorted_bam = params['outPrefix'] + ".sorted.bam"
    
    CMD_1 = CMD_1 % (params['inFastq'], params['outPrefix'], params['index'], params['threads'], params['maxMMap'], params['maxMisMatch'], unsorted_bam)
    CMD_sort_1 = CMD_sort_1 % (unsorted_bam, sorted_bam)
    CMD_sort_2 = CMD_sort_2 % (unsorted_bam, sorted_bam)
    
    print "Start to map to genome:\n\t%s" % (CMD_1, )
    os.system(CMD_1)
    if params['noMut5']:
        print "Start to sort bam:\n\t%s" % (CMD_sort_2, )
        os.system(CMD_sort_2)
    else:
        print "Start to sort bam:\n\t%s" % (CMD_sort_1, )
        os.system(CMD_sort_1)

if __name__ == "__main__":
    main()

