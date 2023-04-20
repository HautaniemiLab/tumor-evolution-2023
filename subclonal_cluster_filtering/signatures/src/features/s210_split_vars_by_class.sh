#!/bin/bash

#Script for splitting variants into separate files by their class (currently
#SBS, DBS, ID). Other types are discarded.

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

#Input and output variables
    IN_DIR="data/interim/filt_set_$NUM_ID/variants_sac_classified/"
    IN_SUF=".tsv"

    SBS_OUT_DIR="data/interim/filt_set_$NUM_ID/variants_sbs/"
    SBS_SUF=".tsv"
    DBS_OUT_DIR="data/interim/filt_set_$NUM_ID/variants_dbs/"
    DBS_SUF=".tsv"
    ID_OUT_DIR="data/interim/filt_set_$NUM_ID/variants_id/"
    ID_SUF=".tsv"

#Make sure output directories exist
    mkdir -p $SBS_OUT_DIR $DBS_OUT_DIR $ID_OUT_DIR

#Clear previous content if it exists
    rm -f $SBS_OUT_DIR/*$SBS_OUT_SUF $DBS_OUT_DIR/*$DBS_OUT_SUF $ID_OUT_DIR/*$ID_OUT_SUF

#Function for filtering variant classes into a file
function filter_variants_by_class() {
    #Input variables
    local FILE=$1
    local IN_SUF=$2
    local OUT_DIR=$3
    local OUT_SUF=$4
    local CLASS=$5

    #Output file
    NEW_FILE=${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF

    #Header
    head -1 $FILE > $NEW_FILE

    #Read variants of specified class
    tail -n +2 $FILE |
        awk '$5 == class' "class=$CLASS" \
        >> $NEW_FILE
}

#Split variants to SBS, DBS, ID
for FILE in $IN_DIR/*$IN_SUF; do
    #Make SBS, DBS and ID files
    filter_variants_by_class $FILE $IN_SUF $SBS_OUT_DIR $SBS_SUF "SBS" &
    filter_variants_by_class $FILE $IN_SUF $DBS_OUT_DIR $DBS_SUF "DBS" &
    filter_variants_by_class $FILE $IN_SUF $ID_OUT_DIR $ID_SUF "ID" &

    wait
done

