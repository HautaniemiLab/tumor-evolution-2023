#!/bin/bash

#Script to type DBS variants i.e. assign mutation type and strand to same
#classification as in Alexandrov et al. ? (bioRxiv 2018)

#Options
    #Defaults
    NUM_ID="0"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )     shift
                                NUM_ID=$1
        esac
        shift
    done

#Name variables
    TYPING_SCRIPT="src/features/s3101_dbs_typing.pl"
    
    IN_DIR="data/interim/filt_set_$NUM_ID/variants_dbs/"
    IN_SUF=".tsv"
    
    OUT_DIR="data/features/filt_set_$NUM_ID/variants_dbs/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Patient-wise function to get mutation type and strand
function patient_dbs_types() {
    #Variables
    local FILE=$1
    local IN_SUF=$2
    local OUT_SUF=$3
    local PATIENT=$(echo $FILE | sed s/$IN_SUF$// | sed 's,'^"$IN_DIR"',,')
    local NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Header
    head -1 $FILE |
        awk '{$5 = $5 "\ttype\tstrand"; print}' OFS='\t' \
        > $NEW_FILE

    #Get type and strand as new columns as final file
    tail -n +2 $FILE |
        $TYPING_SCRIPT \
        >> $NEW_FILE
}

#Loop over files to determine strand and mutation typ
for FILE in $IN_DIR/*$IN_SUF; do
    #DBS with VCF/BAM filtering
    patient_dbs_types $FILE $IN_SUF $OUT_SUF
done

wait

