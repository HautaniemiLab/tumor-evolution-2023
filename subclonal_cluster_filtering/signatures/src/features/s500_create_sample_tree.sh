#!/bin/bash

#Performs hierarchical bottom-up clustering on samples based on number of shared
#variants. Assigns each variant to a tree branch or -2 / -1 / 0, the latter
#where
# -2: variant's depth is too low in all samples
# -1: maximum variant ALT read count among samples is too low
#  0: variant cannot be assigned to a tree branch (ambiguous)

#Trunk can be further split by median VAF to subsections.

#Options
    #Defaults
    NUM_ID="0"
    REGEX_NORMALS=""
    MIN_RC_PATIENT=""
    MIN_RC_SAMPLE=""
    MIN_DEPTH=""
    ROOT_CUTOFFS=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )         shift
                                    NUM_ID=$1 ;;
            -r | --regex-normals )  shift
                                    REGEX_NORMALS="-r $1" ;;
            -min-rc-patient | --minimum-readcount-patient ) shift
                                                            MIN_RC_PATIENT="-min-rc-patient $1" ;;
            -min-rc-sample | --minimum-readcount-sample )   shift
                                                            MIN_RC_SAMPLE="-min-rc-sample $1" ;;
            -min-depth | --minimum-depth )                  shift
                                                            MIN_DEPTH="-min-depth $1" ;;
            -cuts | --root-cutoffs )                        shift
                                                            ROOT_CUTOFFS="-cuts $1"
        esac
        shift
    done

#Name variables
    CLUSTERING_SCRIPT="src/features/s501_sample_tree.R"

    IN_DIR="data/interim/filt_set_$NUM_ID/variants_sac/"
    IN_SUF=".tsv"

    OUT_DIR="data/interim/filt_set_$NUM_ID/sample_tree/"
    OUT_SUF=".RData"

#Make sure target folder exist
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Create sample trees for each patient
for FILE in $IN_DIR/*$SUFFIX; do
    #New file name
    NEW_FILE=${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF

    #Perform clustering
    $CLUSTERING_SCRIPT -i $FILE -o $NEW_FILE $REGEX_NORMALS \
        $MIN_RC_PATIENT $MIN_RC_SAMPLE $MIN_DEPTH $ROOT_CUTOFFS &
done

wait

