#!/bin/bash

#Counts mutation subtype frequencies for SBS (and DBS, ID if available), i.e.
#computes the spectra, for sample tree branches and supersamples which include
#any variant in primary sample tree branches (i.e. 1+).

#Options
    #Defaults
    NUM_ID="0"
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id ) shift
                            NUM_ID=$1 ;;
            -c | --cores )  shift
                            CORES="-c $1"
        esac
        shift
    done

#Name variables
    SPECTRA_SCRIPT="src/features/s511_tree_branch_mutational_spectrum.R"

    IN_DIR_RDATA="data/interim/filt_set_$NUM_ID/sample_tree/"
    IN_SUF_RDATA=".RData"
    IN_SUF_VARS=".tsv"

    OUT_DIR="data/features/filt_set_$NUM_ID/sample_tree/mutational_spectra/"

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
    $SPECTRA_SCRIPT -s $IN_DIR_RDATA -in-st-suf $IN_SUF_RDATA \
        -v $IN_DIR_VARS -in-v-suf $IN_SUF_VARS \
        -m $MUT_CLASS -o $OUT_DIR -out-suf $OUT_SUF $CORES &
done

wait

