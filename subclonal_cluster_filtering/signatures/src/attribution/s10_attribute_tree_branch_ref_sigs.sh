#!/bin/bash

#Performs signature attribution for reference signatures on sample tree super-
#samples and tree branches. This is done for each of the SBS, DBS and ID types.
#
#Attribution done in steps:
# 1. Data set wide forward selection for common signatures
# 2. Supersample backward-forward with patient-wise single tree branch backward-
#    forward followed by supersample backward (w/ rules) for patient/supersample
#    signatures. Only primary branches considered here
# 3. Sample tree branch backward-selection (w/ rules)
#
#Data set and supersample attributions only consider supersamples or branches
#respectively with a minimum number of mutations
#
#In sensitive mode, backward selections are replaced with backward-forward with
#the backward threshold much relaxed in the final supersample backward selection
#step. Affects only multisample patients.
#
#Poisson model optimises generalised Kullback-Leibler instead of cosine
#similarity, and uses BIC (Poisson log-likelihood) to perform model selection.
#Tends to overfit in samples with thousands of mutations, especially with SBS.
#More suitable for WES data.
#
#Data set wide forward selection can be disabled and the set of signatures found
#there or otherwise can be supplemented via options.
#
#Rules refer to SBS rules: forcing SBS1 and SBS5 (default, age-related) and
#connected signatures (e.g. SBS2 and SBS13 present together). Equivalent rules
#are present with COSMIC signatures. Forced signatures can be specified and
#replaced.

#Options
    #Defaults
    NUM_ID="0"
    REFERENCE_SIGNATURES=""
    REGEX_DATASET_SELECTION=""
    NO_DATASET_SELECTION=""
    COMMON_SIGNATURES=""
    FORCED_SIGNATURES=""
    USE_POISSON_MODEL=""
    SENSITIVE_MODE=""
    SENS_BWD_THRES=""
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                         shift
                                                    NUM_ID=$1 ;;
            -r-sel | --regex-dataset-selection )    shift
                                                    REGEX_DATASET_SELECTION="-r-sel $1" ;;
            -ref-sigs | --reference-signatures )    shift
                                                    REFERENCE_SIGNATURES=$1 ;;
            --no-dataset-selection )                NO_DATASET_SELECTION=$1 ;;
            -c-sigs | --common-signatures )         shift
                                                    COMMON_SIGNATURES="-c-sigs $1" ;;
            -f | --forced-signatures )              shift
                                                    FORCED_SIGNATURES="-f $1" ;;
            -p | --use-poisson-model )              USE_POISSON_MODEL=$1 ;;
            --sensitive-mode )                      SENSITIVE_MODE=$1 ;;
            -b | --sensitive-backward-threshold )   shift
                                                    SENS_BWD_THRES="--sensitive-backward-threshold $1" ;;
            -c | --cores )                          shift
                                                    CORES="-c $1"
        esac
        shift
    done

#Name variables
    ATTRIBUTION_SCRIPT="src/attribution/s11_attribute_tree_branch_signatures.R"

    OUT_DIR="results/attributions/filt_set_$NUM_ID/sample_tree/reference/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*.txt $OUT_DIR/*.tsv

#Compute signature attribution for each of the mutation types
for MUT_TYPE in "SBS" "DBS" "ID"; do
    #Name variables
    MUT_TYPE_LOWER=$(echo $MUT_TYPE | tr '/A-Z/' '/a-z/')
    IN_FILE_SUPERSAMPLE="data/features/filt_set_${NUM_ID}/sample_tree/mutational_spectra/supersample_spectra_${MUT_TYPE_LOWER}.tsv"
    IN_FILE_TREE_BRANCH="data/features/filt_set_${NUM_ID}/sample_tree/mutational_spectra/sample_tree_branch_spectra_${MUT_TYPE_LOWER}.tsv"

    #Check if custom reference signature was given
    REF_SIG_PATH="$(echo $REFERENCE_SIGNATURES | perl -pe 's/,/\n/g' |
        perl -e 'while(<STDIN>) {chomp; if (m/^(.*):(.*)$/ and uc $1 eq uc $ARGV[0]) {print $2}}' $MUT_TYPE)"
    if [[ ! -z $REF_SIG_PATH ]]; then
        REF_SIG="-ref-sigs $REF_SIG_PATH"
    else
        REF_SIG=""
    fi

    #Run attribution script
    $ATTRIBUTION_SCRIPT -super $IN_FILE_SUPERSAMPLE -branch $IN_FILE_TREE_BRANCH \
        -m $MUT_TYPE -o $OUT_DIR $REGEX_DATASET_SELECTION \
        $NO_DATASET_SELECTION $COMMON_SIGNATURES $FORCED_SIGNATURES \
        $REF_SIG $USE_POISSON_MODEL $SENSITIVE_MODE $SENS_BWD_THRES $CORES &
done

wait

