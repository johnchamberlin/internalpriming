#!/bin/bash

# pass the following arguments:
#STAR: path to STAR executable
#STARREF: reference index directory
#WHITELIST: 10x genomics whitelist file for V3 chemistry
#CDNA: the cDNA read fastq file
#BCS: the barcodes read fastq file
#OUTDIR: the name of the output directory
# e.g.: /scratch/output
#OUTPREFIX: the name of the full output path with prefix,
# i.e.: /scratch/output/sample1

# FIFO is used to stream the output of seqtk into STARsolo
R1=R1
R2=R2
mkfifo $R1 $R2

# some of the AIBS libraries were mistakenly resequenced to 
# 26bp barcode reads instead of 28bp, and merged into a single fastq
# we need to discard these truncated reads to avoid an error
seqtk mergepe $BCS $CDNA | seqtk seq -L 28 | seqtk dropse \
| tee >(seqtk seq -1 > $R1) | seqtk seq -2 > $R2 &

mkdir $OUTDIR

$STAR --genomeDir $STARREF \
--soloType CB_UMI_Simple --soloCBwhitelist $WHITELIST \
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
--soloBarcodeReadLength 28 \
--readFilesIn $R2 $R1 \
--soloFeatures Gene GeneFull \
--soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts \
--outFileNamePrefix $OUTPREFIX \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--runThreadN 16