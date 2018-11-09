
####################
## Preparation 
####################

cd example_data

make uncompress

## Build index for genome
mkdir index
icSHAPE-pipe starbuild -i fasta/chr22.fa -o index --gtf GTF/chr22.gtf -p 20

## Covert GTF file
icSHAPE-pipe parseGTF -g GTF/chr22.gtf -o GTF/chr22 -s gencode --genome fasta/chr22.fa

## Build index for rRNA
bowtie2-build rRNA/human_rRNA_tRNA_mtRNA.fa rRNA/human_rRNA_tRNA_mtRNA


####################
## Create directory
####################

cd example_data

mkdir -p 1.readCollapse
mkdir -p 2.trim
mkdir -p 3.rem_rRNA
mkdir -p 4.mapGenome
mkdir -p 5.rpkm
mkdir -p 6.sam2tab
mkdir -p 7.calcGenomeSHAPE
mkdir -p 8.transSHAPE
mkdir -p 9.bedGraph
mkdir -p 10.quality_control


####################
## Data process
####################

cd example_data

######## 1. readCollapse

icSHAPE-pipe readcollapse -U raw_data/D1.chr22.fastq.gz -o 1.readCollapse/D1.chr22.fastq
icSHAPE-pipe readcollapse -U raw_data/D2.chr22.fastq.gz -o 1.readCollapse/D2.chr22.fastq
icSHAPE-pipe readcollapse -U raw_data/N1.chr22.fastq.gz -o 1.readCollapse/N1.chr22.fastq
icSHAPE-pipe readcollapse -U raw_data/N2.chr22.fastq.gz -o 1.readCollapse/N2.chr22.fastq

######## 2. Trim reads

icSHAPE-pipe trim -i 1.readCollapse/D1.chr22.fastq -o 2.trim/D1.chr22.fastq -l 13 -a adaptor/TruSeq2-PE.fa -p 20 -m 25
icSHAPE-pipe trim -i 1.readCollapse/D2.chr22.fastq -o 2.trim/D2.chr22.fastq -l 13 -a adaptor/TruSeq2-PE.fa -p 20 -m 25
icSHAPE-pipe trim -i 1.readCollapse/N1.chr22.fastq -o 2.trim/N1.chr22.fastq -l 13 -a adaptor/TruSeq2-PE.fa -p 20 -m 25
icSHAPE-pipe trim -i 1.readCollapse/N2.chr22.fastq -o 2.trim/N2.chr22.fastq -l 13 -a adaptor/TruSeq2-PE.fa -p 20 -m 25

######## 3. Remove rRNA

icSHAPE-pipe cleanFq -i 2.trim/D1.chr22.fastq -o 3.rem_rRNA/D1.chr22.fastq -x rRNA/human_rRNA_tRNA_mtRNA -p 20 --mode End_to_End --sam 3.rem_rRNA/D1.chr22.sam
icSHAPE-pipe cleanFq -i 2.trim/D2.chr22.fastq -o 3.rem_rRNA/D2.chr22.fastq -x rRNA/human_rRNA_tRNA_mtRNA -p 20 --mode End_to_End --sam 3.rem_rRNA/D2.chr22.sam
icSHAPE-pipe cleanFq -i 2.trim/N1.chr22.fastq -o 3.rem_rRNA/N1.chr22.fastq -x rRNA/human_rRNA_tRNA_mtRNA -p 20 --mode End_to_End --sam 3.rem_rRNA/N1.chr22.sam
icSHAPE-pipe cleanFq -i 2.trim/N2.chr22.fastq -o 3.rem_rRNA/N2.chr22.fastq -x rRNA/human_rRNA_tRNA_mtRNA -p 20 --mode End_to_End --sam 3.rem_rRNA/N2.chr22.sam

######## 4. Map to genome

icSHAPE-pipe mapGenome -i 3.rem_rRNA/D1.chr22.fastq -o 4.mapGenome/D1 -x index -p 20 --noMut5
icSHAPE-pipe mapGenome -i 3.rem_rRNA/D2.chr22.fastq -o 4.mapGenome/D2 -x index -p 20 --noMut5
icSHAPE-pipe mapGenome -i 3.rem_rRNA/N1.chr22.fastq -o 4.mapGenome/N1 -x index -p 20 --noMut5
icSHAPE-pipe mapGenome -i 3.rem_rRNA/N2.chr22.fastq -o 4.mapGenome/N2 -x index -p 20 --noMut5

######## 5. Calculate FPKM

icSHAPE-pipe calcFPKM -i 4.mapGenome/D1.sorted.bam -o 5.rpkm/D1 -G GTF/chr22.gtf -p 10
icSHAPE-pipe calcFPKM -i 4.mapGenome/D2.sorted.bam -o 5.rpkm/D2 -G GTF/chr22.gtf -p 10

######## 6. Sam file to tab

icSHAPE-pipe sam2tab -in 4.mapGenome/D1.sorted.bam -out 6.sam2tab/D1.tab
icSHAPE-pipe sam2tab -in 4.mapGenome/D2.sorted.bam -out 6.sam2tab/D2.tab
icSHAPE-pipe sam2tab -in 4.mapGenome/N1.sorted.bam -out 6.sam2tab/N1.tab
icSHAPE-pipe sam2tab -in 4.mapGenome/N2.sorted.bam -out 6.sam2tab/N2.tab

######## 7. Calculate SHAPE score

icSHAPE-pipe calcSHAPE \
    -D 6.sam2tab/D1.tab,6.sam2tab/D2.tab \
    -N 6.sam2tab/N1.tab,6.sam2tab/N2.tab \
    -size index/chrNameLength.txt \
    -ijf index/sjdbList.fromGTF.out.tab \
    -out 7.calcGenomeSHAPE/shape.gTab

######## 8. generate transcriptome-based SHAPE score and RTBD

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i 7.calcGenomeSHAPE/shape.gTab \
    -o 8.transSHAPE/shape.out \
    -g GTF/chr22.genomeCoor.bed \
    -s index/chrNameLength.txt \
    -r 5.rpkm/D1/isoforms.fpkm_tracking,5.rpkm/D2/isoforms.fpkm_tracking

######## 9. Visualization

icSHAPE-pipe genSHAPEToBedGraph -i 7.calcGenomeSHAPE/shape.gTab -t TrtCont -o 9.bedGraph/ -c 200


##########################
##### Replicate control
##########################

######## 1. Reads distribution statistic

icSHAPE-pipe readDistributionStatistic \
    -1 raw_data/D1.chr22.fastq.gz,raw_data/D2.chr22.fastq.gz,raw_data/N1.chr22.fastq.gz,raw_data/N2.chr22.fastq.gz \
    -2 1.readCollapse/D1.chr22.fastq,1.readCollapse/D2.chr22.fastq,1.readCollapse/N1.chr22.fastq,1.readCollapse/N2.chr22.fastq \
    -3 2.trim/D1.chr22.fastq,2.trim/D2.chr22.fastq,2.trim/N1.chr22.fastq,2.trim/N2.chr22.fastq \
    -4 3.rem_rRNA/D1.chr22.fastq,3.rem_rRNA/D2.chr22.fastq,3.rem_rRNA/N1.chr22.fastq,3.rem_rRNA/N2.chr22.fastq \
    -5 4.mapGenome/D1.unsorted.bam,4.mapGenome/D2.unsorted.bam,4.mapGenome/N1.unsorted.bam,4.mapGenome/N2.unsorted.bam \
    --labels DMSO_1,DMSO_2,NAI_1,NAI_2 \
    -o 10.quality_control/read_map_statistics.pdf


######## 2. Reads map statistic

icSHAPE-pipe samStatistics \
    -i 4.mapGenome/D1.unsorted.bam \
    -o 10.quality_control/D1.pdf \
    -t 10.quality_control/D1_report.txt \
    -g GTF/chr22.genomeCoor.bed \
    --fast

icSHAPE-pipe samStatistics \
    -i 4.mapGenome/D2.unsorted.bam \
    -o 10.quality_control/D2.pdf \
    -t 10.quality_control/D2_report.txt \
    -g GTF/chr22.genomeCoor.bed \
    --fast

icSHAPE-pipe samStatistics \
    -i 4.mapGenome/N1.unsorted.bam \
    -o 10.quality_control/N1.pdf \
    -t 10.quality_control/N1_report.txt \
    -g GTF/chr22.genomeCoor.bed \
    --fast

icSHAPE-pipe samStatistics \
    -i 4.mapGenome/N2.unsorted.bam \
    -o 10.quality_control/N2.pdf \
    -t 10.quality_control/N2_report.txt \
    -g GTF/chr22.genomeCoor.bed \
    --fast


######## 3. RT replicates

icSHAPE-pipe countRT \
    -in 6.sam2tab/D1.tab,6.sam2tab/D2.tab,6.sam2tab/N1.tab,6.sam2tab/N2.tab \
    -ijf index/sjdbList.fromGTF.out.tab \
    -size index/chrNameLength.txt \
    -out 10.quality_control/RT_count.txt


icSHAPE-pipe plotGenomeRTRepCor -i 10.quality_control/RT_count.txt -o 10.quality_control/DMSO_rt_cor.pdf --col1 4 --col2 6
icSHAPE-pipe plotGenomeRTRepCor -i 10.quality_control/RT_count.txt -o 10.quality_control/NAI_rt_cor.pdf --col1 8 --col2 10

######## 4. SHAPE replicates

icSHAPE-pipe calcSHAPE \
    -D 6.sam2tab/D1.tab \
    -N 6.sam2tab/N1.tab \
    -size index/chrNameLength.txt \
    -ijf index/sjdbList.fromGTF.out.tab \
    -out 10.quality_control/shape_rep1.gTab

icSHAPE-pipe calcSHAPE \
    -D 6.sam2tab/D2.tab \
    -N 6.sam2tab/N2.tab \
    -size index/chrNameLength.txt \
    -ijf index/sjdbList.fromGTF.out.tab \
    -out 10.quality_control/shape_rep2.gTab

icSHAPE-pipe combine_gTab_SHAPE 10.quality_control/shape_rep1.gTab 10.quality_control/shape_rep2.gTab 10.quality_control/shape_combine.txt
icSHAPE-pipe plotGenomeSHAPERepCor -i 10.quality_control/shape_combine.txt -o 10.quality_control/shape_cor.pdf

######## 5. Evaluate icSHAPE with known 18S structure

icSHAPE-pipe sam2tab -in 3.rem_rRNA/D1.chr22.sam -out 10.quality_control/D1.rRNA.tab
icSHAPE-pipe sam2tab -in 3.rem_rRNA/D2.chr22.sam -out 10.quality_control/D2.rRNA.tab
icSHAPE-pipe sam2tab -in 3.rem_rRNA/N1.chr22.sam -out 10.quality_control/N1.rRNA.tab
icSHAPE-pipe sam2tab -in 3.rem_rRNA/N2.chr22.sam -out 10.quality_control/N2.rRNA.tab

icSHAPE-pipe calcSHAPE \
    -D 10.quality_control/D1.rRNA.tab,10.quality_control/D2.rRNA.tab \
    -N 10.quality_control/N1.rRNA.tab,10.quality_control/N2.rRNA.tab \
    -size rRNA/human_rRNA_tRNA_mtRNA.len \
    -out 10.quality_control/human_rRNA_tRNA_mtRNA.gTab

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i 10.quality_control/human_rRNA_tRNA_mtRNA.gTab \
    -s rRNA/human_rRNA_tRNA_mtRNA.len \
    -o 10.quality_control/human_rRNA_tRNA_mtRNA.out

icSHAPE-pipe evaluateSHAPE \
    -i 10.quality_control/human_rRNA_tRNA_mtRNA.out \
    -s rRNA/human_18S.dot \
    -o 10.quality_control/human_18S.pdf



