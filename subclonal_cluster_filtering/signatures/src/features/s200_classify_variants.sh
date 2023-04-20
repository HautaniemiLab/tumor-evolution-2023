#!/bin/bash

#Adds a column of mutation class for variants:
# - SBS (Single Base Substitution)
# - ID (InDel)
# - DBS (Double Base Substitution)
# - MBS (Multiple Base Substitution)
# - other (none of the above; not following correct format)

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
    IN_DIR="data/interim/filt_set_$NUM_ID/variants_sac/"
    IN_SUF=".tsv"

    OUT_DIR="data/interim/filt_set_$NUM_ID/variants_sac_classified/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Loop over files
for FILE in $IN_DIR/*$IN_SUF; do
    #New file name
    NEW_FILE=${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF

    #Header
    head -1 $FILE |
        awk '{$4 = $4 "\tclass"; print}' OFS='\t' \
        > $NEW_FILE

    #Classify variants based on REF and ALT lengths
    tail -n +2 $FILE |
        awk '{
            if ((length($3) == 1 && length($4) > 1) || (length($3) > 1 && length($4) == 1) ) {
                $4 = $4 "\tID"
            } else if (length($3) == length($4) && length($4) == 1) {
                $4 = $4 "\tSBS"
            } else if (length($3) == length($4) && length($4) == 2) {
                $4 = $4 "\tDBS"
            } else if (length($3) == length($4) && length($4) > 2) {
                $4 = $4 "\tMBS"
            } else {
                $4 = $4 "\tother"
            }
            print
        }' OFS='\t' \
        >> $NEW_FILE

done

