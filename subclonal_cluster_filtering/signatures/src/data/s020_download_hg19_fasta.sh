#!/bin/bash

#Downloads hg19 reference genome fasta, unpacks and combines to a single file
#and finally indexes it.

#File name variables
    IN_FILE="hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"

    OUT_DIR="data/raw/reference_genomes/"
    TEMP_CHR_DIR="${OUT_DIR}.tmp.hg19_chromosomes/"

    OUT_FILE="${OUT_DIR}hg19.fa"

#Make the directories exist
    mkdir -p $OUT_DIR
    mkdir -p $TEMP_CHR_DIR

#Download
    wget -nv $IN_FILE -O $OUT_FILE

#Uncopmress
    tar -zxf $OUT_FILE -C $TEMP_CHR_DIR

#Create combined fasta
    cat ${TEMP_CHR_DIR}chr*.fa > $OUT_FILE

#Index the output
    samtools faidx $OUT_FILE

#Remove the compressed genome and the temporary fastas
    rm -r $TEMP_CHR_DIR

