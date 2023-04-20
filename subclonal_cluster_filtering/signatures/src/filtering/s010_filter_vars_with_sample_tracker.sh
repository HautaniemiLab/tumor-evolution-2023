#!/bin/bash

#Filter variants from a patient such that no variants discovered only in excluded
#samples are kept. Uses the filtering script and tracker table from the variants
#pipeline. Modifies the specified sample list to the format required by the
#filtering script.
#The sample list can be a proper metadata table which can be used via filtering
#arguments (R/dplyr style consecutive) to choose samples to retain.
#The filtering is turned off, i.e. skipped, by default. Instead, input is simply
#copied to the output.

#Options
    #Defaults
    NUM_ID="0"
    SKIP_SAMPLE_FILTER="1"
    TRACKER_TABLE="data/raw/sample_tracker.csv"
    SAMPLE_LIST=""
    SAMPLE_COL_NAME="sample"
    SAMPLE_FILTERS_CONCATENATED=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )             shift
                                        NUM_ID=$1 ;;
            -t | --tracker-table )      shift
                                        TRACKER_TABLE=$1 ;;
            -s | --sample-list )        shift
                                        SAMPLE_LIST=$1 ;;
            --sample-column-name )      shift
                                        SAMPLE_COL_NAME="$1" ;;
            --sample-filters-concatenated ) shift
                                            SAMPLE_FILTERS_CONCATENATED="$1 $SAMPLE_FILTERS_CONCATENATED" ;;
            -sf | --sample-filter )         shift
                                            SAMPLE_FILTERS_CONCATENATED="$SAMPLE_FILTERS_CONCATENATED -sf $1" ;;
            --sample-tracker-filter )   SKIP_SAMPLE_FILTER="0"
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
    FILTERING_SCRIPT="data/raw/filter_variants_by_sample.pl"

    IN_DIR="data/interim/filt_set_$NUM_ID/variants_tsv_annot/"
    IN_SUF=".tsv"

    OUT_DIR="data/interim/filt_set_$NUM_ID/variants_filt_by_sample/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Process sample list to single line comma-separated string (samples come from
#specific column name
#If no valid sample list is given (or default does not exist), skip filtering
    if [[ -f $SAMPLE_LIST ]]; then
        #If sample column name is not in the table header, table is assumed
        #headerless, and sample names assumed to be on the first column. In such
        #a case, all filtering, besides the table's listed samples, is ignored
        SAMPLES_KEPT=$(
            grep -vP "^#|^>|^\-|^$" $SAMPLE_LIST |
                R --slave -e '
                    suppressMessages(require(dplyr))

                    args = commandArgs(trailingOnly=TRUE)

                    sample_col_name = args[1]
                    filters_concatenated = args[2]

                    sample_list = read.table("stdin", header=T, sep="\t", as.is=T)

                    if (sample_col_name %in% colnames(sample_list)) {
                        for (filter_string in strsplit(filters_concatenated, " -sf ")[[1]][-1]) {
                            sample_list = eval(parse(
                                text = paste("sample_list %>%", filter_string)
                            ))
                        }

                        cat(paste(sample_list[,sample_col_name], collapse=","))
                    } else {
                        cat(
                            paste(
                                c(colnames(sample_list)[1], sample_list[,1]),
                                collapse = ","
                            )
                        )
                    }
                ' --args $SAMPLE_COL_NAME "$SAMPLE_FILTERS_CONCATENATED"
        )
    else
        echo "No valid sample list \"$SAMPLE_LIST\"; skipping sample filter"
        SKIP_SAMPLE_FILTER="1"
    fi

#Function for deriving new file name
function new_file_name() {
    local IN_FILE=$1
    local IN_SUF=$2
    local OUT_SUF=$3

    echo ${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF
}

#Function for processing a patient CSV
function filter_patient_csv() {
    #Input variable
    local IN_FILE=$1

    #Output file
    local OUT_FILE=$(new_file_name $IN_FILE $IN_SUF $OUT_SUF)

    #Operate the filtering script
    $FILTERING_SCRIPT -v $IN_FILE -s $SAMPLES_KEPT -t $TRACKER_TABLE \
        -o $OUT_FILE

    #Remove empty results (patient filtered)
    if [ $(wc -l $OUT_FILE | sed 's/\s.*$//') -eq 1 ]; then
        rm $OUT_FILE
    fi
}

#Function for simply copying the input to target
function copy_csv() {
    #Input variable
    local IN_FILE=$1

    #Output file
    local OUT_FILE=$(new_file_name $IN_FILE $IN_SUF $OUT_SUF)

    cp $IN_FILE $OUT_FILE
}

#Process all patient files
for FILE in $IN_DIR/*$IN_SUF; do
    #Simply copy the file if skipping sample filtering
    if [[ $SKIP_SAMPLE_FILTER == 0 ]]; then
        filter_patient_csv $FILE &
    else
        copy_csv $FILE
    fi
done

wait

