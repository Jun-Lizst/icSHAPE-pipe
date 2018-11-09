#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import commands
import gzip
import version

Usage = """
readCollapse - Remove duplicate reads from fastq
================================================
\x1b[1mUSAGE:\x1b[0m
  %s -1 fastq_PE_reads_1 -2 fastq_PE_reads_2 -U fastq_SE_reads
\x1b[1mHELP:\x1b[0m
  -U     Single ends read
  -1     Paired ends read 1
  -2     Paired ends read 2

 More options
  -o     Single ends read output 
  -p     PE read output 1
  -q     PE read output 2

  -f     Unique fasta after collapse
  -l     <1-10> barcode length (default: guess)
         Longer barcode, less memory, more time

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

dirname = os.path.dirname(os.path.abspath(__file__))
SplitBarcode = "python " + os.path.join(dirname, 'splitByBarcode.py')
collapseFq = "python " + os.path.join(dirname, 'collapseSingleFq.py')


def init():
    import getopt
    
    Params = { 'outFasta': None, 'bcLen': None }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hU:1:2:o:p:q:f:l:')
    
    for op, value in opts:
        if op == '-h':
            print Usage
            exit(-1)
        # Basic Parameters
        elif op == '-U':
            Params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFile'] = os.path.abspath(value)
        
        elif op == '-1':
            Params['inFile_1'] = os.path.abspath(value)
        elif op == '-2':
            Params['inFile_2'] = os.path.abspath(value)
        elif op == '-p':
            Params['outFile_1'] = os.path.abspath(value)
        elif op == '-q':
            Params['outFile_2'] = os.path.abspath(value)
        
        elif op == '-l':
            Params['bcLen'] = int(value)
            assert 1 <= Params['bcLen'] <= 10
        
        elif op == '-f':
            Params['outFasta'] = os.path.abspath(value)
        
        else:
            print >>sys.stderr, "parameter Error: unrecognized parameter: "+op
            print Usage
            sys.exit(-1)
    
    if 'inFile' in Params:
        Params['mode'] = 'single'
        print >>sys.stderr, "Single-end mode..."
        if 'outFile' not in Params:
            print >>sys.stderr, "Error: \n\tSigle-end mode: specify -U, -o\n\tPair-end mode: specify -1, -2, -p and -q"
            exit(-1)
        outFile, outDir, fileSuffix = fileparse(Params['outFile'])
        Params['outDir'] = outDir
    else:
        Params['mode'] = 'pair'
        print >>sys.stderr, "Pair-end mode..."
        if ('inFile_1' not in Params) or ('inFile_2' not in Params) or ('outFile_1' not in Params) or ('outFile_2' not in Params):
            print >>sys.stderr, "Error: \n\tSigle-end mode: specify -U, -o\n\tPair-end mode: specify -1, -2, -p and -q"
            exit(-1)
        outFile, outDir, fileSuffix = fileparse(Params['outFile_1'])
        Params['outDir'] = outDir
    
    return Params

def main():
    Params = init()
    
    tmpDir = Params['outDir'] + "/tmp_"+str(os.getpid())
    if not os.path.exists(tmpDir):
        print >>sys.stderr, commands.getoutput("mkdir "+tmpDir)
    
    if Params['mode'] == 'single':
        print >>sys.stderr, "Collapsing file %s...\n\t" % (Params['inFile'], );
        if os.path.exists(Params['outFile']):
            print >>sys.stderr, "Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile'], );
            print >>sys.stderr, commands.getoutput("rm "+Params['outFile'])
        inFile = Params['inFile']
        outFile = Params['outFile']
    else:
        if os.path.exists(Params['outFile_1']):
            print >>sys.stderr, "Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile_1'], );
            print >>sys.stderr, commands.getoutput("rm "+Params['outFile_1'])
        if os.path.exists(Params['outFile_2']):
            print >>sys.stderr, "Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile_2'], );
            print >>sys.stderr, commands.getoutput("rm "+Params['outFile_2'])
        inFile = tmpDir + "/input.fastq"
        outFile = tmpDir + "/tmpOut.fastq"
        mergePairEndReads(Params['inFile_1'], Params['inFile_2'], inFile)
    
    if Params['outFasta'] and os.path.exists(Params['outFasta']):
        print >>sys.stderr, "Warning! %s exisits, will be overwritten. \n\t" % (Params['outFasta'], );
        print >>sys.stderr, commands.getoutput("rm "+Params['outFasta'])
    
    file2Collapse = [ inFile ]
    if Params['bcLen'] != None:
        cPos = 5
        cLen = Params['bcLen']
    else:
        cPos, cLen = estimateSplit(inFile)
    
    #print cPos, cLen
    if cLen>0:
        print >>sys.stderr, "File %s too large, will be splitted..." % ( inFile,  )
        file2Collapse.pop()
        tmpOutDir = tmpDir
        CMD = SplitBarcode+" -i %s -o %s -p %s -l %s --mode new" % (inFile, tmpOutDir, cPos, cLen)
        code = os.system(CMD)
        if code != 0:
            print >>sys.stderr, "Error! splitting file failed. quiting..."
            exit(-1)
        else:
            file2Collapse = getFqFiles( tmpOutDir )
            #print file2Collapse
    
    totalReads = 0
    uniqReads = 0
    outputFastas = []
    for file in file2Collapse:
        collapseResults = ""
        if Params['outFasta'] != None:
            CMD = collapseFq+" -i %s -o %s --mode append --fasta %s" % (file, outFile, Params['outFasta'])
        else:
            CMD = collapseFq+" -i %s -o %s --mode append" % (file, outFile)
        #print CMD
        code = os.system(CMD)
        #volReads, volUniqReads, volUniqRatio, volFasta = _parseCollapseOutput ( "......." );
        #totalReads += volReads
        #uniqReads += volUniqReads
        #outputFastas.append( volFasta )
    
    #uniqRatio = ".2f" % (1.0*uniqReads/totalReads)
    if Params['mode'] == 'pair':
        splitPairEndReads ( outFile, Params['outFile_1'], Params['outFile_2'] )
        print >>sys.stderr, commands.getoutput("rm -f "+inFile)
        print >>sys.stderr, commands.getoutput("rm -f "+outFile)
        print >>sys.stderr,  "Collapsing file %s and %s finished." % (Params['inFile_1'], Params['inFile_2'])
    else:
        print >>sys.stderr,  "Collapsing file %s finished." % (Params['inFile'], )
    
    #print >>sys.stderr, "Read collapse successful! Total count: %s, unique count: %s, unique ratio: %s." % (totalReads, uniqReads, uniqRatio)
    os.system("rm -r "+tmpDir)

def mergePairEndReads(readFile1, readFile2, peFile):
    ## should test whether they are of the same length
    print >>sys.stderr, "merge two PE fastq files..."
    inFile1 = readFile1
    inFile2 = readFile2
    if readFile1.endswith('.gz'):
        inFile1 = readFile1 + ".tmp.fq"
        print >>sys.stderr, commands.getoutput("gzip -d -c %s > %s" % (readFile1, inFile1))
    elif readFile2.endswith('.gz'):
        inFile2 = readFile2 + ".tmp.fq"
        print >>sys.stderr, commands.getoutput("gzip -d -c %s > %s" % (readFile2, inFile2))
    
    print >>sys.stderr, commands.getoutput("paste %s %s > %s" % (inFile1, inFile2, peFile))

def estimateSplit(inputFile):
    import math
    
    fileSize = os.path.getsize(inputFile)
    
    pos = 1
    length = int( math.log(fileSize/1000000000.0, 4.0) )
    if inputFile.endswith('.gz'):
        IN = gzip.open(inputFile, 'r')
    else:
        IN = open(inputFile, 'r')
    line = IN.readline; line = IN.readline()
    IN.close()
    readLen = len(line)
    if readLen <= 10:
        length = 0
    else:
        pos = int( readLen/4.0 + 10 )
        while pos + length > readLen:
            pos -= length
            if pos < 1:
                pos = 1
                break
    
    return pos, length

def splitPairEndReads(peFile, readFile1, readFile2):
    print >>sys.stderr, "split into two PE fastq files..."
    PE = gzip.open(peFile) if peFile.endswith('.gz') else open(peFile)
    R1 = gzip.open(readFile1, 'w') if readFile1.endswith('.gz') else open(readFile1, 'w')
    R2 = gzip.open(readFile2, 'w') if readFile2.endswith('.gz') else open(readFile2, 'w')
    for line in PE:
        r1, r2 = line.split('\t')
        R1.writelines(r1)
        R2.writelines(r2)
    PE.close()
    R1.close()
    R2.close()

def getFqFiles(Dir):
    files = os.listdir(Dir)
    fastqFiles = [ Dir.rstrip('/')+"/"+file for file in files if file.endswith(".fastq") or file.endswith(".fq") ]
    return fastqFiles

def fileparse(fillFilePath):
    fileDir = '/'.join(fillFilePath.split('/')[:-1])
    fileName = fillFilePath.split('/')[-1]
    pureFileName = '.'.join(fileName.split('.')[:-1])
    fileSuffix = fileName.split('.')[-1]
    return pureFileName, fileDir, fileSuffix

if __name__ == '__main__':
    main()

