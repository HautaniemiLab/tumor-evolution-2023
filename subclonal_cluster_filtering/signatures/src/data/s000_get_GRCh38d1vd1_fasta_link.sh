#!/bin/bash

#Creates soft link of the reference genome (GRCh38d1vd1) from resources. Also
#makes sure to index the file.

#File name variables
    IN_FILE="/path/to/GRCh38.d1.vd1.fa"
    OUT_DIR="data/raw/reference_genomes/"
    OUT_FILE="${OUT_DIR}GRCh38.d1.vd1.fa"

#Make sure the output directory exists
    mkdir -p $OUT_DIR

#Create link
    ln -sf $IN_FILE $OUT_FILE

#Index the fasta file
    samtools faidx $OUT_FILE

