#!/bin/bash

#Draws sample trees in a VAF plot with variants ordered by tree branches.

#Options
    #Defaults
    NUM_ID="0"
    SORT_BY_VAF=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )         shift
                                    NUM_ID=$1 ;;
            -s | --sort-by-vaf )    SORT_BY_VAF=$1
        esac
        shift
    done

#Name variables
    PLOT_SCRIPT="src/visualisation/s31_vaf_tree_tracks_plot.R"

    IN_DIR="data/interim/filt_set_$NUM_ID/sample_tree/"
    IN_SUF=".RData"

    OUT_DIR="results/plots/filt_set_$NUM_ID/sample_tree/"
    OUT_SUF=".pdf"

#Make sure target folder exist
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Draw sample trees for each patient
for FILE in $IN_DIR/*$SUFFIX; do
    #New file name
    NEW_FILE=${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF

    #Draw plot
    $PLOT_SCRIPT -i $FILE -o $NEW_FILE $SORT_BY_VAF &
done

wait

