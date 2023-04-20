#!/bin/bash

#Counts mutation subtype frequencies for SBS (and DBS, ID if available), i.e.
#computes the spectra, for single samples (at least one ALT read) and
#supersamples which include any variant with at least 1 ALT read in tumours.

#Options
    #Defaults
    NUM_ID="0"
    REGEX_NORMALS=""
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )         shift
                                    NUM_ID=$1 ;;
            -r | --regex-normals )  shift
                                    REGEX_NORMALS="-r $1" ;;
            -c | --cores )          shift
                                    CORES="-c $1"
        esac
        shift
    done

#Name variables
    SPECTRA_SCRIPT="src/features/s401_mutational_spectrum.R"

    IN_SUFFIX=".tsv"

    OUT_DIR="data/features/filt_set_$NUM_ID/mutational_spectra/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*.tsv

#Compute the spectra for each mutation class
for MUT_CLASS in "SBS" "DBS" "ID"; do
    #Name variables
    MUT_CLASS_LOWER=$(echo $MUT_CLASS | tr '/A-Z/' '/a-z/')
    IN_DIR="data/features/filt_set_$NUM_ID/variants_${MUT_CLASS_LOWER}/"
    OUT_SUF="_${MUT_CLASS_LOWER}.tsv"

    #Create the mutational spectra
    $SPECTRA_SCRIPT -i $IN_DIR -in-suf $IN_SUFFIX -m $MUT_CLASS \
        -o $OUT_DIR -out-suf $OUT_SUF $REGEX_NORMALS $CORES &
done

wait

