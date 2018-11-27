#-*- coding:utf-8 -*-


def init(genomeCoorBedFile, seqFn='', showAttr=True, rem_tVersion=False, rem_gVersion=False):
    import ParseTrans
    return ParseTrans.ParseTransClass(genomeCoorBedFile, seqFileName=seqFn, showAttr=showAttr, remove_tid_version=rem_tVersion, remove_gid_version=rem_gVersion)



def initGTF(AnnotationGTF, genomeFile='', source='Gencode', showAttr=True, rem_tVersion=False, rem_gVersion=False, verbose=False):
    """
    Build a parser from GTF or GFF3 file
    source: Gencode or NCBI
    """
    
    if source == 'Gencode':
        import GENCODE_Genome
        handle = GENCODE_Genome.GENCODE_Genome_Class(AnnotationGTF)
        Parser = __build_parser(handle, genomeFn=genomeFile, source='Gencode', showAttr=showAttr, rem_tVersion=rem_tVersion, rem_gVersion=rem_gVersion, verbose=verbose)
    elif source == 'NCBI':
        import NCBI_Genome
        handle = NCBI_Genome.NCBI_Genome_Class(AnnotationGTF)
        Parser = __build_parser(handle, genomeFn=genomeFile, source='NCBI', showAttr=showAttr, rem_tVersion=rem_tVersion, rem_gVersion=rem_gVersion, verbose=verbose)
    else:
        print >>sys.stderr, "Error: source must be Gencode or NCBI"
    
    return Parser


## Edit on 2018-11-10
def __build_parser(GenomeHandler, genomeFn='', source='Gencode', showAttr=True, rem_tVersion=False, rem_gVersion=False, verbose=False):
    import random, os, commands, ParseTrans
    rID = random.randint(10000, 99999)
    tmp_genomeCoor_file = "/tmp/tmp_%s_genomeCoor.bed" % (rID, )
    
    if source == 'Gencode':
        GenomeHandler.write_genomeCoor_bed(tmp_genomeCoor_file, onlyChr=False, pureTransID=rem_tVersion, pureGeneID=rem_gVersion, verbose=verbose)
    elif source == 'NCBI':
        GenomeHandler.write_genomeCoor_bed(tmp_genomeCoor_file, onlyChr=False, pureTransID=rem_tVersion, verbose=verbose)
    else:
        print >>sys.stderr, "Error!"
        return None
    
    tmp_transcriptome_file = ""
    if genomeFn:
        tmp_transcriptome_file = "/tmp/tmp_%s_transcriptome.fasta" % (rID, )
        GenomeHandler.writeTranscriptome(genomeFn, tmp_transcriptome_file, pureTransID=rem_tVersion, onlyChr=False, verbose=verbose)
    
    Parser = ParseTrans.ParseTransClass(tmp_genomeCoor_file, seqFileName=tmp_transcriptome_file, showAttr=showAttr, remove_tid_version=rem_tVersion, remove_gid_version=rem_gVersion)
    
    os.remove(tmp_genomeCoor_file)
    if tmp_transcriptome_file:
        os.remove(tmp_transcriptome_file)
    
    return Parser

