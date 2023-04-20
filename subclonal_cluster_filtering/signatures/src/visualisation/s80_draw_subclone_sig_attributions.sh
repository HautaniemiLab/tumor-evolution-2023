#!/bin/bash

#Draws exposure and triple barplots for subclonal supersample and cluster
#signature attributions. Exposure plots for clusters include indicator of the
#cluster and their sizes.

#Options
    #Defaults
    NUM_ID="0"
    REFERENCE_SIGNATURES=""
    ZERO_INDEXED_CLUSTERS=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                         shift
                                                    NUM_ID=$1 ;;
            -ref-sigs | --reference-signatures )    shift
                                                    REFERENCE_SIGNATURES=$1 ;;
            -z | --zero-indexed-clusters )          ZERO_INDEXED_CLUSTERS="-z"
        esac
        shift
    done

#Name variables
    PLOT_SCRIPT="src/visualisation/s81_subclone_sig_attr_plotting.R"

    IN_DIR="data/interim/filt_set_$NUM_ID/subclones/"
    IN_SUF=".RData"

    OUT_DIR="results/plots/filt_set_$NUM_ID/subclones/attributions/reference/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -rf $OUT_DIR/*

#Compute signature attribution for each of the mutation types
for MUT_TYPE in "SBS" "DBS" "ID"; do
    #Name variables
    MUT_TYPE_LOWER=$(echo $MUT_TYPE | tr '/A-Z/' '/a-z/')
    IN_FILE_CATAL_SUPERSAMPLE="data/features/filt_set_${NUM_ID}/subclones/mutational_spectra/supersample_spectra_${MUT_TYPE_LOWER}.tsv"
    IN_FILE_CATAL_CLUSTERS="data/features/filt_set_${NUM_ID}/subclones/mutational_spectra/subclone_spectra_${MUT_TYPE_LOWER}.tsv"
    IN_FILE_EXPS_SUPERSAMPLE="results/attributions/filt_set_${NUM_ID}/subclones/reference/contributions_${MUT_TYPE_LOWER}_supersample.tsv"
    IN_FILE_EXPS_CLUSTERS="results/attributions/filt_set_${NUM_ID}/subclones/reference/contributions_${MUT_TYPE_LOWER}_subclone.tsv"

    #Check if custom reference signature was given
    REF_SIG_PATH="$(echo $REFERENCE_SIGNATURES | perl -pe 's/,/\n/g' |
        perl -e 'while(<STDIN>) {chomp; if (m/^(.*):(.*)$/ and uc $1 eq uc $ARGV[0]) {print $2}}' $MUT_TYPE)"
    if [[ ! -z $REF_SIG_PATH ]]; then
        REF_SIG="-ref-sigs $REF_SIG_PATH"
    else
        REF_SIG=""
    fi

#TODO Handle this better in general
    if [[ -f $IN_FILE_EXPS_SUPERSAMPLE && -f $IN_FILE_EXPS_CLUSTERS ]]; then
        #Run attribution script
        $PLOT_SCRIPT -super-catal $IN_FILE_CATAL_SUPERSAMPLE -subclone-catal $IN_FILE_CATAL_CLUSTERS \
            -super-exp $IN_FILE_EXPS_SUPERSAMPLE -subclone-exp $IN_FILE_EXPS_CLUSTERS \
            -m $MUT_TYPE -o $OUT_DIR $REF_SIG $ZERO_INDEXED_CLUSTERS &
    else
        echo "Skipping $MUT_TYPE mutations as contribution tables are missing!"
    fi
done

wait

