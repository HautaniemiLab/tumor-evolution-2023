#!/bin/bash

#Downloads hs37d5 reference genome fasta, unpacks it and indexes it.

#File name variables
    IN_FILE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    OUT_DIR="data/raw/reference_genomes/"
    OUT_FILE="${OUT_DIR}hs37d5.fa"

#Make sure the output directory exists
    mkdir -p $OUT_DIR

#Download
    wget -nv $IN_FILE -O ${OUT_FILE}.gz

#Uncopmress
    gunzip ${OUT_FILE}.gz

#Index the resulting fasta file
    samtools faidx $OUT_FILE

