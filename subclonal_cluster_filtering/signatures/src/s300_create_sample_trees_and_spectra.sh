#!/bin/bash

#Mini shell workflow for creating crude sample trees by hierarchically
#clustering variants among samples of a patient and then computing mutational
#spectra for the sample trees as well as their combined supersample.
#Plots sample tree VAF plots and mutational profiles.

#Save full command
FULL_COMMAND="$BASH_SOURCE $@"

#Options
    #Defaults
    NUM_ID="0"
    REGEX_NORMALS=""
    MIN_RC_PATIENT="2"
    MIN_RC_SAMPLE="1"
    MIN_DEPTH="10"
    ROOT_CUTOFFS=""
    SORT_BY_VAF=""
    PLOT_SEPARATE="-p"
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                                 shift
                                                            NUM_ID=$1 ;;
            -r | --regex-normals )                          shift
                                                            REGEX_NORMALS="-r $1" ;;
            -min-rc-patient | --minimum-readcount-patient ) shift
                                                            MIN_RC_PATIENT=$1 ;;
            -min-rc-sample | --minimum-readcount-sample )   shift
                                                            MIN_RC_SAMPLE=$1 ;;
            -min-depth | --minimum-depth )                  shift
                                                            MIN_DEPTH=$1 ;;
            -cuts | --root-cutoffs )                        shift
                                                            ROOT_CUTOFFS="-cuts $1" ;;
            -s | --sort-by-vaf )                            SORT_BY_VAF=$1 ;;
            --do-not-plot-separate )                        PLOT_SEPARATE="" ;;
            -c | --cores )                                  shift
                                                            CORES="-c $1"
        esac
        shift
    done

#Write version, uncommitted changes and full command to troubleshooting file
    src/troubleshoot/log_command.sh --command "$FULL_COMMAND" \
        --log-file "data/troubleshoot/filt_set_$NUM_ID/commands/$(basename $BASH_SOURCE)"

#Print used options
    echo "Analysing variants with sample filter set \"$NUM_ID\""
    echo "- The following filtering thresholds are used:"
    echo "  - Minimum depth in all samples $MIN_DEPTH"
    echo "  - ALT read count in patient (at least one sample) $MIN_RC_PATIENT and in sample $MIN_RC_SAMPLE"
    if [[ $REGEX_NORMALS != "" ]]; then
        echo "- Regex string \"$(echo $REGEX_NORMALS | perl -pe 's/^-r ?//')\" is used to match normal samples"
    fi
    if [[ $ROOT_CUTOFFS != "" ]]; then
        echo "- Splitting tree root with $(echo $ROOT_CUTOFFS | perl -pe 's/^-cuts ?//') cutoffs"
    fi
    if [[ $SORT_BY_VAF == "-s" ]]; then
        echo "- Variants within VAF tree branches will be sorted by median VAF"
    fi
    if [[ $PLOT_SEPARATE == "-p" ]]; then
        echo "- Tree branch specific mutational profiles are plotted to separate files by patient"
    fi
    echo ""

#Sample tree
    #Creation
    echo "Creating sample trees"
    src/features/s500_create_sample_tree.sh -n $NUM_ID $REGEX_NORMALS \
        -min-rc-patient $MIN_RC_PATIENT -min-rc-sample $MIN_RC_SAMPLE -min-depth $MIN_DEPTH \
        $ROOT_CUTOFFS

    #Plotting VAF heatmap
    echo "Plotting sample tree VAF heatmaps"
    src/visualisation/s30_draw_vaf_tree.sh -n $NUM_ID $SORT_BY_VAF

#Mutational spectra
    #Creation
    echo "Creating mutational spectra of SBS, DBS and ID mutations"
    src/features/s510_create_tree_branch_spectra.sh -n $NUM_ID $CORES

    #Plotting
    echo "Plotting mutational spectra"
    src/visualisation/s40_draw_tree_branch_spectra.sh -n $NUM_ID $PLOT_SEPARATE

