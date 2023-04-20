#!/bin/bash

#Downloads GENCODE annotations and unpacks it.

#File name variables
    IN_FILE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
    OUT_DIR="data/raw/genome_annotations/"
    OUT_FILE="${OUT_DIR}GRCh38.GENCODE.v29.gtf.gz"

#Download
    wget -nv $IN_FILE -O $OUT_FILE

#Uncopmress
    gunzip $OUT_FILE

