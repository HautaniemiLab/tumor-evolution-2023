#!/bin/bash

#Counts mutation subtype frequencies for SBS (and DBS, ID if available), i.e.
#computes the spectra, for subclonal clusters and supersamples which include
#variants from all clusters.

#Options
    #Defaults
    NUM_ID="0"
    CLUSTERS_DIR="data/features/subclones/"
    SUBCLONE_INFO="data/raw/subclonal_clusters.csv"
    SUBCLONE_INFO=""
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                     shift
                                                NUM_ID=$1 ;;
            -i | --input-clusters-dir )         shift
                                                CLUSTERS_DIR=$1 ;;
            -s-info | --subclone-info )         shift
                                                SUBCLONE_INFO="-s-info $1" ;;
            -c | --cores )                      shift
                                                CORES="-c $1"
        esac
        shift
    done

#Name variables
    SPECTRA_SCRIPT="src/features/s621_subclone_mutational_spectrum.R"

    IN_DIR_CLUSTERS="$CLUSTERS_DIR/cluster_variants/"
    IN_SUF_CLUSTERS=".tsv"
    IN_SUF_VARS=".tsv"

    OUT_DIR="data/features/filt_set_$NUM_ID/subclones/mutational_spectra/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*.tsv

#Compute the spectra for each mutation class
for MUT_CLASS in "SBS" "DBS" "ID"; do
    #Name variables
    MUT_CLASS_LOWER=$(echo $MUT_CLASS | tr '/A-Z/' '/a-z/')
    IN_DIR_VARS="data/features/filt_set_$NUM_ID/variants_${MUT_CLASS_LOWER}/"
    OUT_SUF="_${MUT_CLASS_LOWER}.tsv"

    #Create the mutational spectra
    $SPECTRA_SCRIPT -s $IN_DIR_CLUSTERS -in-s-suf $IN_SUF_CLUSTERS \
        -v $IN_DIR_VARS -in-v-suf $IN_SUF_VARS -m $MUT_CLASS \
        -o $OUT_DIR -out-suf $OUT_SUF $SUBCLONE_INFO $CORES &
done

wait

