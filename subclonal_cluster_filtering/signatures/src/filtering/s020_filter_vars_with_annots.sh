#!/bin/bash

#Multipurpose filtering script processes the variant call csv files to filter
#columns, eeping only variant identifying CHR, POS, REF, ALT as well as SAC for
#each sample specified to be retained (via a sample list and associated filters).
#Also clears quotes and removes NA lines if any.
#If SAC columns are missing, they will be computed from F1R2 and F2R1.

#The sample list can be a proper metadata table which can be used via filtering
#arguments (R/dplyr style consecutive) to choose samples to retain.

#Additionally filters at variant level; variants must satisfy:
# - #ALT reads > 0, ignores normal
# - Options for minimum + and - strand ALT read count and minimum VAF (in at
#   least one sample), ignores normal
# - Options for additional free R/dplyr style consecutive variant filters using
#   any annotations in the variant tables
# - Optional maximum VAF table with tumour sample name (low purity sample to
#   filter likely germline mutations)

#The max VAF filter is applied before sample filters; others are applied after.

#TODO: perhaps should be combined with s000_filter_samples.sh to reduce
#      redundancy, in particular with sample filtering code
#      - Problem is support for old data sets requiring the separate tracker
#        functionality

#Options
    #Defaults
    NUM_ID="0"
    SAMPLE_LIST=""
    SAMPLE_COL_NAME="sample"
    VAF_THRES="0.05"
    FWD_THRES="1"
    REV_THRES="1"
    SAMPLE_FILTERS_CONCATENATED=""
    VARIANT_FILTERS_CONCATENATED=""
    REGEX_NORMALS=""
    VAF_FILTER_TABLE=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )             shift
                                        NUM_ID=$1 ;;
            -s | --sample-list )        shift
                                        SAMPLE_LIST=$1 ;;
            --vaf )                     shift
                                        VAF_THRES=$1 ;;
            --fwd )                     shift
                                        FWD_THRES=$1 ;;
            --rev )                     shift
                                        REV_THRES=$1 ;;
            --sample-column-name )      shift
                                        SAMPLE_COL_NAME="$1" ;;
            --sample-filters-concatenated )     shift
                                                SAMPLE_FILTERS_CONCATENATED="$1 $SAMPLE_FILTERS_CONCATENATED" ;;
            -sf | --sample-filter )             shift
                                                SAMPLE_FILTERS_CONCATENATED="$SAMPLE_FILTERS_CONCATENATED -sf $1" ;;
            --variant-filters-concatenated )    shift
                                                VARIANT_FILTERS_CONCATENATED="$1 $VARIANT_FILTERS_CONCATENATED" ;;
            -vf | --variant-filter )            shift
                                                VARIANT_FILTERS_CONCATENATED="$VARIANT_FILTERS_CONCATENATED -vf $1" ;;
            -r | --regex-normals )      shift
                                        REGEX_NORMALS=$1 ;;
            --vaf-filter-table )        shift
                                        VAF_FILTER_TABLE="--vaf-filter-table $1"
        esac
        shift
    done

    #Default sample lists; one derived from num-id or without num-id in order should they exist
    if [[ $SAMPLE_LIST == "" ]]; then
        if [[ -f "data/sample_filter_list_$NUM_ID" ]]; then
            SAMPLE_LIST="data/sample_filter_list_$NUM_ID"
            echo "Using default sample list \"$SAMPLE_LIST\"!"
        elif [[ -f "data/sample_filter_list" ]]; then
            SAMPLE_LIST="data/sample_filter_list"
            echo "Using default sample list \"$SAMPLE_LIST\"!"
        else
            echo "Neither default sample list paths exist: \"data/sample_filter_list_$NUM_ID\" or \"data/sample_filter_list\"!"
        fi
    fi

#Input and output variables
    FILTERING_SCRIPT="src/filtering/s021_filter_samples_vars_and_columns.R"

    IN_DIR="data/interim/filt_set_$NUM_ID/variants_filt_by_sample/"
    IN_SUF=".tsv"

    OUT_DIR="data/interim/filt_set_$NUM_ID/variants_sac/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Process all patient files
for FILE in $IN_DIR/*$IN_SUF; do
    #Name variables
    PATIENT=`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] /x)' $FILE $IN_SUF`
    NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Apply filtering
    $FILTERING_SCRIPT -i $FILE -s $SAMPLE_LIST -o $NEW_FILE \
        --sample-column-name $SAMPLE_COL_NAME \
        --sample-filters-concatenated "$SAMPLE_FILTERS_CONCATENATED" \
        --variant-filters-concatenated "$VARIANT_FILTERS_CONCATENATED" \
        -r $REGEX_NORMALS --vaf $VAF_THRES --fwd $FWD_THRES --rev $REV_THRES \
        $VAF_FILTER_TABLE --patient $PATIENT &
done

wait

