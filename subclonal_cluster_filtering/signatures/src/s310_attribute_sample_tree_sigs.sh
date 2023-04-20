#!/bin/bash

#Mini shell workflow for attributing reference signatures to mutational
#spectra/profiles for SBS, DBS and ID mutations of both supersamples (full
#trees) and single sample tree branches.
#Plots drawn include signature contributions and fits.

#Save full command
FULL_COMMAND="$BASH_SOURCE $@"

#Options
    #Defaults
    NUM_ID="0"
    REFERENCE_SIGNATURES=""
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
            -ref-sigs | --reference-signatures )    shift
                                                    REFERENCE_SIGNATURES="-ref-sigs $1" ;;
            -r-sel | --regex-dataset-selection )    shift
                                                    REGEX_DATASET_SEL="-r-sel $1" ;;
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

#Write version, uncommitted changes and full command to troubleshooting file
    src/troubleshoot/log_command.sh --command "$FULL_COMMAND" \
        --log-file "data/troubleshoot/filt_set_$NUM_ID/commands/$(basename $BASH_SOURCE)"

#Print used options
    echo "Analysing variants with sample filter set \"$NUM_ID\""
    if [[ $REGEX_DATASET_SELECTION != "" ]]; then
        echo "- Regex string \"$(echo $REGEX_DATASET_SELECTION | perl -pe 's/^-r ?//')\" is used to match patients not considered in selecting common signatures in the data set"
    fi
    if [[ $REFERENCE_SIGNATURES != "" ]]; then
        echo "- The following custom reference signature paths have been specified:"
        for LINE in $(echo $REFERENCE_SIGNATURES | perl -pe 's/^-ref-sigs //; s/,/\n/g'); do
            echo -n "  - "
            echo -e $(echo $LINE | perl -ne 's/(.*):/:\\t/; print uc($1) . $_')
        done
    fi
    if [[ $NO_DATASET_SELECTION != "" ]]; then
        echo "- Not performing data set level forward selection to obtain common signatures"
    fi
    if [[ $COMMON_SIGNATURES != "" ]]; then
        echo "- The following signatures will be added to common signatures: $(echo $COMMON_SIGNATURES | perl -pe 's/^-c-sigs ?//')"
    fi
    if [[ $FORCED_SIGNATURES != "" ]]; then
        echo "- The following custom forced signatures have been set: $(echo $FORCED_SIGNATURES | perl -pe 's/^-f ?//')"
    else
        echo "- Using default forced signatures (SBS1,SBS5 for sigProfiler, Signatures 1,5 for COSMIC)"
    fi
    if [[ $SENSITIVE_MODE != "" ]]; then
        echo "- Using more sensitive model selection in attribution with multisample patients"

        if [[ $SENS_BWD_THRES != "" ]]; then
            echo "  - Using custom single-to-supersample backward threshold: $(echo $SENS_BWD_THRES | perl -pe 's/^-b ?//')"
        fi
    fi
    if [[ $USE_POISSON_MODEL != "" ]]; then
        echo "- Attribution will use Poisson-based model, optimising generalised Kullback-Leibler. Model selection with BIC (Poisson)"
    else
        echo "- Attribution will use default model, optimising cosine similarity. Model selection with accuracy i.e. cosine similarity"
    fi
    echo ""

#Signature attribution
    echo "Attributing signatures"
    src/attribution/s10_attribute_tree_branch_ref_sigs.sh -n $NUM_ID $REGEX_DATASET_SELECTION \
        $REFERENCE_SIGNATURES $NO_DATASET_SELECTION $COMMON_SIGNATURES $FORCED_SIGNATURES \
        $USE_POISSON_MODEL $SENSITIVE_MODE $SENS_BWD_THRES $CORES

#Plot signature exposures and fits
    echo "Plotting signature exposures and fits"
    src/visualisation/s50_draw_tree_branch_sig_attributions.sh -n $NUM_ID $REFERENCE_SIGNATURES

