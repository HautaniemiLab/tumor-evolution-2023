#!/bin/bash

#Mini shell workflow for filtering variants data, classifying variants and
#computing the mutational profiles/spectra for SBS, DBS and ID mutations, with
#one for supersamples as well as single samples (variant exists if there is 1
#ALT read). Steps include:
# - Converting VCF to CSV (tab-separated) with requested annotations (useful
#   for filtering), with optional blacklist or whitelist filtering; can be
#   skipped if using CSV input. Input files can be skipped with a grep
#   filename (basename) regex
# - Filtering by variant discovery in samples (turned off by default; option
#   for old somatic variant versions; deprecated)
# - Filtering samples using a sample filtering table (metadata), with options
#   to use metadata fields to choose/remove samples (R/dplyr style consecutive)
# - Filtering variants by read counts, VAF and optionally variant annotations
# - Classifying and splitting variants (SBS/DBS/ID)
# - Subtyping variants (spectra subtypes)
# - Computing profiles (subtype frequency)
# - Plotting the profiles

#Save full command
FULL_COMMAND="$BASH_SOURCE $@"

#Options
    #Defaults
    NUM_ID="0"
    VAF_THRES="0.05"
    FWD_THRES="1"
    REV_THRES="1"
    IN_DIR="data/raw/variants_vcf/"
    INFO_QUERY_STRING=""
    FORMAT_QUERY_STRING="\t%AF\t%F1R2\t%F2R1"
    TRACKER_TABLE=""
    SAMPLE_LIST=""
    SAMPLE_COL_NAME=""
    SAMPLE_FILTERS_CONCATENATED=""
    VARIANT_FILTERS_CONCATENATED=""
    SKIP_SAMPLE_FILTER=""
    PLOT_SEPARATE="-p"
    REGEX_NORMALS=""
    VAF_FILTER_TABLE=""
    FILENAME_REGEX=""
    BLACKLIST=""
    INVERT_BLACKLIST=""
    FILTER_MBS_AS_SBS=""
    DONT_KEEP_NESTED_SBS=""
    CSV=""
    CORES=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )             shift
                                        NUM_ID=$1 ;;
            --vaf )                     shift
                                        VAF_THRES=$1 ;;
            --fwd )                     shift
                                        FWD_THRES=$1 ;;
            --rev )                     shift
                                        REV_THRES=$1 ;;
            -i | --input-directory  )   shift
                                        IN_DIR=$1 ;;
            --do-not-plot-separate )    PLOT_SEPARATE="" ;;
            --info-query )              shift
                                        INFO_QUERY_STRING=$1 ;;
            --format-query )            shift
                                        FORMAT_QUERY_STRING=$1 ;;
            -t | --tracker-table )      shift
                                        TRACKER_TABLE="-t $1"
                                        SAMPLE_TRACKER_FILTER="--sample-tracker-filter" ;;
            -s | --sample-list )        shift
                                        SAMPLE_LIST="-s $1" ;;
            --sample-column-name )      shift
                                        SAMPLE_COL_NAME="--sample-column-name $1" ;;
            -sf | --sample-filter )     shift
                                        SAMPLE_FILTERS_CONCATENATED="$SAMPLE_FILTERS_CONCATENATED -sf $1" ;;
            -vf | --variant-filter )    shift
                                        VARIANT_FILTERS_CONCATENATED="$VARIANT_FILTERS_CONCATENATED -vf $1" ;;
            --sample-tracker-filter )   SAMPLE_TRACKER_FILTER=$1 ;;
            -r | --regex-normals )      shift
                                        REGEX_NORMALS="-r $1" ;;
            --vaf-filter-table )        shift
                                        VAF_FILTER_TABLE="--vaf-filter-table $1" ;;
            --filename-regex )          shift
                                        FILENAME_REGEX="$1" ;;
            --blacklist )               shift
                                        BLACKLIST="--blacklist $1" ;;
            --invert-blacklist )        INVERT_BLACKLIST=$1 ;;
            --filter-mbs-as-sbs )       FILTER_MBS_AS_SBS=$1 ;;
            --dont-keep-nested-sbs )    DONT_KEEP_NESTED_SBS=$1 ;;
            --csv )                     CSV=$1 ;;
            -c | --cores )              shift
                                        CORES="-c $1"
        esac
        shift
    done

#Write version, uncommitted changes and full command to troubleshooting file
    src/troubleshoot/log_command.sh --command "$FULL_COMMAND" \
        --log-file "data/troubleshoot/filt_set_$NUM_ID/commands/$(basename $BASH_SOURCE)"

#Print used options
    #NUM_ID and input directory
    echo "Analysing variants with sample filter set \"$NUM_ID\", from input directory \"$IN_DIR\""

    #File regex filter
    if [[ $FILENAME_REGEX != "" ]]; then
        echo "- Filtering input files using \"grep $(echo $FILENAME_REGEX | perl -pe 's/^--filename-regex //')\""
    fi

    #CSV or VCF INFOs and FORMATs
    if [[ $CSV == "--csv" ]]; then
        echo "- Using CSV input"
    else
        echo "- Using VCF input with extra INFO query \"$INFO_QUERY_STRING\" and FORMAT query \"$FORMAT_QUERY_STRING\""

        #Blacklist
        if [[ $BLACKLIST != "" ]]; then
            if [[ $INVERT_BLACKLIST == "--invert-blacklist" ]]; then
                echo "- Using $(echo $BLACKLIST | perl -pe 's/^--blacklist //') as variant whitelist in VCF-to-CSV conversion"
            else
                echo "- Using $(echo $BLACKLIST | perl -pe 's/^--blacklist //') as variant blacklist in VCF-to-CSV conversion"
            fi
            if [[ $FILTER_MBS_AS_SBS == "--filter-mbs-as-sbs" ]]; then
                echo "  - Filtering MBSs also considers component SBSs in VCF-to-CSV conversion"
                if [[ $DONT_KEEP_NESTED_SBS == "--dont-keep-nested-sbs" ]]; then
                    echo "    - Nested SBS surviving the MBS filter won't be retained"
                fi
            fi
        fi
    fi

    #Sample tracker table
    if [[ $SAMPLE_TRACKER_FILTER == "--sample-tracker-filter" ]]; then
        echo "- Variant filtering by sample discovery has been activated;"
        if [[ $TRACKER_TABLE != "" ]]; then
            echo "  - using custom tracker table \"$(echo $TRACKER_TABLE | perl -pe 's/^-t //')\" to filter variants"
        else
            echo "  - using default tracker table to filter variants"
        fi
    fi

    #Sample list / metadata table
    if [[ $SAMPLE_LIST != "" ]]; then
        echo "- Using custom sample list \"$(echo $SAMPLE_LIST | perl -pe 's/^-s //')\" to filter samples"
    else
        echo "- Using default sample list \"data/sample_filter_list_$NUM_ID\" or alternatively \"data/sample_filter_list\" to filter samples"
    fi

    #Sample filters
    if [[ $SAMPLE_FILTERS_CONCATENATED != "" ]]; then
        echo -n "- The following R/dplyr filters are used to filter samples using the sample (metadata) list, if it is valid:"
        echo "$SAMPLE_FILTERS_CONCATENATED" | sed 's/\ -sf\ /\n  - /g'
    fi

    #Variant filters
    echo "- The following standard variant filtering thresholds are used:"
    echo "  - VAF=$VAF_THRES, minimum forward and reverse reads $FWD_THRES,$REV_THRES"
    if [[ $VARIANT_FILTERS_CONCATENATED != "" ]]; then
        echo -n "- The following R/dplyr filters are used to filter variants:"
        echo "$VARIANT_FILTERS_CONCATENATED" | sed 's/\ -vf\ /\n  - /g'
    fi
    if [[ $VAF_FILTER_TABLE != "" ]]; then
        echo "- Using maximum VAF threshold table \"$(echo $VAF_FILTER_TABLE | perl -pe 's/^--vaf-filter-table //')\" to filter likely germline mutations"
    fi

    #Other
    if [[ $REGEX_NORMALS != "" ]]; then
        echo "- Regex string \"$(echo $REGEX_NORMALS | perl -pe 's/^-r //')\" is used to match normal samples"
    fi
    if [[ $PLOT_SEPARATE == "-p" ]]; then
        echo "- Sample-specific mutational profiles are plotted to separate files by patient"
    fi
    echo ""

#Variant filtering
    #Convert VCF to CSV
    if [[ $CSV != "--csv" ]]; then
        echo "Converting VCFs to CSVs with annotations"
    fi
    src/filtering/s000_convert_vcf_to_csv_with_annots.sh -n $NUM_ID -i $IN_DIR \
        $(if [[ $FILENAME_REGEX != "" ]]; then echo -n "--filename-regex"; fi) "$FILENAME_REGEX" \
        $(if [[ $INFO_QUERY_STRING != "" ]]; then echo -n "--info-query"; fi) "$INFO_QUERY_STRING" \
        $(if [[ $FORMAT_QUERY_STRING != "" ]]; then echo -n "--format-query"; fi) "$FORMAT_QUERY_STRING" \
        $BLACKLIST $INVERT_BLACKLIST $FILTER_MBS_AS_SBS $DONT_KEEP_NESTED_SBS $CSV

    #Filter by discovery in samples
    if [[ $SAMPLE_TRACKER_FILTER == "--sample-tracker-filter" ]]; then
        echo "Filtering variants by discovery in samples"
    fi
    src/filtering/s010_filter_vars_with_sample_tracker.sh -n $NUM_ID \
        $SAMPLE_LIST $SAMPLE_COL_NAME $SAMPLE_TRACKER_FILTER $TRACKER_TABLE \
        --sample-filters-concatenated "$SAMPLE_FILTERS_CONCATENATED" \

    #Filter by FS, SOR, ALT reads and VAF
    echo "Filtering samples then variants by ALT reads, VAF and variant annotations"
    src/filtering/s020_filter_vars_with_annots.sh -n $NUM_ID \
        --vaf $VAF_THRES --fwd $FWD_THRES --rev $REV_THRES \
        $SAMPLE_LIST $SAMPLE_COL_NAME $VAF_FILTER_TABLE $REGEX_NORMALS \
        --sample-filters-concatenated "$SAMPLE_FILTERS_CONCATENATED" \
        --variant-filters-concatenated "$VARIANT_FILTERS_CONCATENATED" \

#Classify variants and split them
    #Classify
    echo "Classifying variants"
    src/features/s200_classify_variants.sh -n $NUM_ID

    #Split
    echo "Splitting variants by class"
    src/features/s210_split_vars_by_class.sh -n $NUM_ID

#Subtype variants
    #SBS contexts
    echo "Subtyping SBS variants with trinucleotide contexts"
    src/features/s3000_sbs_get_contexts.sh -n $NUM_ID

    #SBS transcriptional strand
    echo "Getting SBS transcriptional strand"
    src/features/s3010_sbs_get_ts_status.sh -n $NUM_ID

    #DBS
    echo "Subtyping DBS variants"
    src/features/s3100_dbs_typing.sh -n $NUM_ID

    #ID
    echo "Subtyping ID variants"
    src/features/s3200_id_typing.sh -n $NUM_ID

#Create mutational spectra
    echo "Creating mutational spectra of SBS, DBS and ID mutations"
    src/features/s400_create_mutational_spectra.sh -n $NUM_ID $REGEX_NORMALS $CORES

#Plot mutational spectra
    echo "Plotting mutational spectra"
    src/visualisation/s00_draw_mutational_spectra.sh -n $NUM_ID $PLOT_SEPARATE $REGEX_NORMALS

