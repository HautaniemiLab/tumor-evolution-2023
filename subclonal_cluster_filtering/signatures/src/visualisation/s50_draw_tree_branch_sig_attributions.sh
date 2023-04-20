#!/bin/bash

#Draws exposure and triple barplots for sample tree supersample and branch
#signature attributions. Exposure plots also include dendrograms.

#Options
    #Defaults
    NUM_ID="0"
    REFERENCE_SIGNATURES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                         shift
                                                    NUM_ID=$1 ;;
            -ref-sigs | --reference-signatures )    shift
                                                    REFERENCE_SIGNATURES=$1
        esac
        shift
    done

#Name variables
    PLOT_SCRIPT="src/visualisation/s51_tree_branch_sig_attr_plotting.R"

    IN_DIR="data/interim/filt_set_$NUM_ID/sample_tree/"
    IN_SUF=".RData"

    OUT_DIR="results/plots/filt_set_$NUM_ID/sample_tree/attributions/reference/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -rf $OUT_DIR/*

#Compute signature attribution for each of the mutation types
for MUT_TYPE in "SBS" "DBS" "ID"; do
    #Name variables
    MUT_TYPE_LOWER=$(echo $MUT_TYPE | tr '/A-Z/' '/a-z/')
    IN_FILE_CATAL_SUPERSAMPLE="data/features/filt_set_${NUM_ID}/sample_tree/mutational_spectra/supersample_spectra_${MUT_TYPE_LOWER}.tsv"
    IN_FILE_CATAL_TREE_BRANCH="data/features/filt_set_${NUM_ID}/sample_tree/mutational_spectra/sample_tree_branch_spectra_${MUT_TYPE_LOWER}.tsv"
    IN_FILE_EXPS_SUPERSAMPLE="results/attributions/filt_set_${NUM_ID}/sample_tree/reference/contributions_${MUT_TYPE_LOWER}_supersample.tsv"
    IN_FILE_EXPS_TREE_BRANCH="results/attributions/filt_set_${NUM_ID}/sample_tree/reference/contributions_${MUT_TYPE_LOWER}_sample_tree_branch.tsv"

    #Check if custom reference signature was given
    REF_SIG_PATH="$(echo $REFERENCE_SIGNATURES | perl -pe 's/,/\n/g' |
        perl -e 'while(<STDIN>) {chomp; if (m/^(.*):(.*)$/ and uc $1 eq uc $ARGV[0]) {print $2}}' $MUT_TYPE)"
    if [[ ! -z $REF_SIG_PATH ]]; then
        REF_SIG="-ref-sigs $REF_SIG_PATH"
    else
        REF_SIG=""
    fi

#TODO Handle this better in general
    if [[ -f $IN_FILE_EXPS_SUPERSAMPLE && -f $IN_FILE_EXPS_TREE_BRANCH ]]; then
        #Run attribution script
        $PLOT_SCRIPT -in-st-dir $IN_DIR -in-st-suf $IN_SUF \
            -super-catal $IN_FILE_CATAL_SUPERSAMPLE -branch-catal $IN_FILE_CATAL_TREE_BRANCH \
            -super-exp $IN_FILE_EXPS_SUPERSAMPLE -branch-exp $IN_FILE_EXPS_TREE_BRANCH \
            -m $MUT_TYPE -o $OUT_DIR $REF_SIG &
    else
        echo "Skipping $MUT_TYPE mutations as contribution tables are missing!"
    fi
done

wait

