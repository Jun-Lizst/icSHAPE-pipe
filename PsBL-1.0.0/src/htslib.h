
#ifndef HTSLIB_H
#define HTSLIB_H

#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>
#include "pan_type.h"

namespace pan{

using std::string;
using std::to_string;


/*
SRR3194440.15842008 32  ENST00000527779.1   495 1   5S10M100N10M6S  *   0   0   TTGTGAAGACATTGGTGTTGAACCTGAAAAT JJIHIJJJGJJJJJJHHHIIHIIJJJIJJJJ AS:i:0  XS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:31 YT:Z:UU
*/


string getBamQName(bam1_t *record);
uint16_t getBamFlag(bam1_t *record);
string getBamRef(bam1_t *record, bam_hdr_t *hdr);
int32_t getBamRefPos(bam1_t *record);
int32_t getBamMapQuanlity(bam1_t *record);
string getBamCigar(bam1_t *record);
string getBamMateRef(bam1_t *record, bam_hdr_t *hdr);
int32_t getBamMateRefPos(bam1_t *record);
/*** TLEN skip ***/
string getBamSeq(bam1_t *record);
string getBamQuanlity(bam1_t *record);
string getBamTag(bam1_t *record);


bool getBamHead(bam_hdr_t *hdr, MapStringT<uLONG> &chr_len);

};
#endif


