#!/bin/bash

#Draws sample trees in a VAF plot with variants ordered by tree branches and VAF
#with subclonal cluster tracks.

#Options
    #Defaults
    NUM_ID="0"
    CLUSTERS_DIR="data/features/subclones/"
    SORT_BY_VAF=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )             shift
                                        NUM_ID=$1 ;;
            -i | --input-clusters-dir ) shift
                                        CLUSTERS_DIR=$1 ;;
            -s | --sort-by-vaf )        SORT_BY_VAF=$1
        esac
        shift
    done

#Input and output variables
    PLOT_SCRIPT="src/visualisation/s31_vaf_tree_tracks_plot.R"

    IN_DIR_RDATA="data/interim/filt_set_$NUM_ID/sample_tree/"
    IN_SUF_RDATA=".RData"
    IN_DIR_CLUSTERS="$CLUSTERS_DIR/cluster_variants/"
    IN_SUF_CLUSTERS=".tsv"

    OUT_DIR="results/plots/filt_set_$NUM_ID/subclones/sample_tree/"
    OUT_SUF=".pdf"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Draw plot for each patient
for FILE in $IN_DIR_CLUSTERS/*$IN_SUF_CLUSTERS; do
    #Name variables
    PATIENT=`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF_CLUSTERS`
    RDATA_FILE=${IN_DIR_RDATA}${PATIENT}$IN_SUF_RDATA
    NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Skip visualisation if sample tree data is missing i.e. have not been defined
    if [[ -f $RDATA_FILE ]]; then
        $PLOT_SCRIPT -i $RDATA_FILE -o $NEW_FILE -t $FILE $SORT_BY_VAF &
    fi
done

wait

