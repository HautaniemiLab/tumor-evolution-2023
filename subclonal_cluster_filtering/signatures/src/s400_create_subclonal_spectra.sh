#!/bin/bash

#Mini shell workflow for computing mutational spectra for subclonal mutation
#clusters as well as combined supersamples and clusters from select subclones.
#Plots sample tree VAF plots with subclonal mutations annotated and mutational
#profiles.

#Save full command
FULL_COMMAND="$BASH_SOURCE $@"

#Options
    #Defaults
    NUM_ID="0"
    SUBCLONES_DIR="data/raw/subclones/"
    SUBCLONES_SUFFIX=""
    SUBCLONES_INFO="data/raw/subclones.csv"
    SORT_BY_VAF="-s"
    PLOT_SEPARATE="-p"
    SUBCLONES_INFO=""
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )                 shift
                                            NUM_ID=$1 ;;
            -i | --input-directory )        shift
                                            SUBCLONES_DIR=$1 ;;
            -s-info | --subclone-info )     shift
                                            SUBCLONES_INFO="-s-info $1" ;;
            -s-suffix | --subclone-suffix ) shift
                                            SUBCLONES_SUFFIX="-s-suffix $1" ;;
            --do-not-sort-by-vaf)           SORT_BY_VAF="" ;;
            --do-not-plot-separate )        PLOT_SEPARATE="" ;;
            -c | --cores )                  shift
                                            CORES="-c $1"
        esac
        shift
    done

#Write version, uncommitted changes and full command to troubleshooting file
    src/troubleshoot/log_command.sh --command "$FULL_COMMAND" \
        --log-file "data/troubleshoot/filt_set_$NUM_ID/commands/$(basename $BASH_SOURCE)"

#Print used options
    echo "Analysing variants from sample filter set \"$NUM_ID\""
    echo "- Using subclones from directory \"$SUBCLONES_DIR\""
    if [[ $SUBCLONES_SUFFIX != "" ]]; then
        echo "  - Using custom \"$(echo $SUBCLONES_SUFFIX | perl -pe 's/^-s-suffix //')\" subclone file suffix"
    fi
    if [[ $SUBCLONES_INFO == "" ]]; then
        echo "- Not filtering subclonal clusters"
    else
        echo "- Using specified subclones for supersamples and combined clusters from table \"$SUBCLONES_INFO\""
    fi
    if [[ $SORT_BY_VAF == "-s" ]]; then
        echo "- Variants within VAF tree branches will be sorted by median VAF"
    fi
    if [[ $PLOT_SEPARATE == "-p" ]]; then
        echo "- Subclone specific mutational profiles are plotted to separate files by patient"
    fi
    echo ""

#Generate subclonal data
    #Extract subclonal data
    echo "Extracting subclonal variants"
    src/features/s600_prepare_subclone_tracks.sh -i $SUBCLONES_DIR $SUBCLONES_SUFFIX

    CLUSTER_VARIANTS_DIR="data/features/$(basename $SUBCLONES_DIR)/"

    #Plotting VAF heatmap
    echo "Plotting sample tree VAF heatmaps"
    src/visualisation/s60_draw_vaf_tree_with_subclones.sh -n $NUM_ID \
        -i $CLUSTER_VARIANTS_DIR $SORT_BY_VAF

#Mutational spectra
    #Join variants to subclones and create spectra
    echo "Creating mutational spectra of SBS, DBS and ID mutations"
    src/features/s620_create_subclone_spectra.sh -n $NUM_ID \
        -i $CLUSTER_VARIANTS_DIR $SUBCLONES_INFO $CORES

    #Plot mutational spectra
    echo "Plotting mutational spectra"
    src/visualisation/s70_draw_subclone_spectra.sh -n $NUM_ID $PLOT_SEPARATE

