#!/bin/bash

#Convert VCF variants to headered TSV with given annotations (INFO and FORMAT
#filtered).
#Uses arguments to control INFO and FORMAT annotations to be output.
#Input files can be filtered using a filename (basename) regex filter via
#grep. This regex can be inverted by prefixing with '-v '.
#It also comes with an option to filter variants via a VCF blacklist. It can
#be turned into a whitelist via an option.
#The conversion can be skipped if the input is CSV via argument --csv.

#Options
    #Defaults
    NUM_ID="0"
    IN_DIR="data/raw/variants_vcf/"
    INFO_QUERY_STRING=""
    FORMAT_QUERY_STRING="\t%AF\t%F1R2\t%F2R1"
    FILENAME_REGEX=""
    BLACKLIST=""
    INVERT_BLACKLIST="0"
    FILTER_MBS_AS_SBS="0"
    CSV="0"

    IN_SUF=".vcf"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )             shift
                                        NUM_ID=$1 ;;
            -i | --input-directory )    shift
                                        IN_DIR=$1 ;;
            --info-query )              shift
                                        INFO_QUERY_STRING=$1 ;;
            --format-query )            shift
                                        FORMAT_QUERY_STRING=$1 ;;
            --filename-regex )          shift
                                        FILENAME_REGEX=$1 ;;
            --blacklist )               shift
                                        BLACKLIST=$1 ;;
            --invert-blacklist )        INVERT_BLACKLIST="1" ;;
            --filter-mbs-as-sbs )       FILTER_MBS_AS_SBS="1" ;;
            --dont-keep-nested-sbs )    DONT_KEEP_NESTED_SBS=$1 ;;
            --csv )                     CSV="1"
                                        IN_SUF=".csv"
        esac
        shift
    done

#Output variables
    OUT_DIR="data/interim/filt_set_$NUM_ID/variants_tsv_annot/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*

#Function for deriving new file name
function new_file_name() {
    local IN_FILE=$1
    local IN_SUF=$2
    local OUT_SUF=$3

    echo ${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] (\.gz)? $ /x)' $IN_FILE $IN_SUF`$OUT_SUF
}

#Function for filtering VCF from variants in a blacklist
function filter_vcf_with_blacklist() {
    #Input variables
    local IN_FILE=$1
    local OUT_FILE=$2

    #Make sure input is VCF.gz with indexing
    bcftools view $IN_FILE -Oz -o "$OUT_FILE.temp.vcf.gz"
    tabix -p vcf "$OUT_FILE.temp.vcf.gz"

    #Should MBS be filtered using their SBS components?
    if [[ $FILTER_MBS_AS_SBS == 1 ]]; then
#TODO: MBS as SBS filtering and inversion does not work properly (need to rethink logic)
#      - The original annotated MBS must survive the isec phase (they should bypass first isec)...
#      - Apply the MBS w/ SBS filter
#      - Apply isec on the MBS
#      - It should work this way without inversion when the blacklist DOES contain MBS
#      - Inversion logic: should an MBS be retained if any one nested SBS or full MBS is in the blacklist?
#        (currently requiring ALL)...
#      - Also note that filtered MBS would output smaller MBS / SBS (depending on argument)...
#      - Maybe isec => split MBS | others => isec (nested SBS only) | - => join MBS | - => merge
        trap "rm -rf $OUT_FILE.temp.vcf.gz $OUT_FILE.temp.vcf.gz.tbi $OUT_FILE.temp.mbs_split.vcf.gz $OUT_FILE.temp.mbs_split.vcf.gz.tbi $OUT_FILE.temp.isec.vcf.gz" RETURN

        #Add SBS of MBSs
        src/filtering/s001_split_mbs_to_sbs_with_id.sh -i $OUT_FILE.temp.vcf.gz -o $OUT_FILE.temp.mbs_split.vcf.gz

        #Apply blacklist filter
        bcftools isec -c none -n=$(expr 1 + $INVERT_BLACKLIST) -w1 $OUT_FILE.temp.mbs_split.vcf.gz $BLACKLIST \
            -Oz -o $OUT_FILE.temp.isec.vcf.gz

        #Figure out MBS/SBS filtering, output appropriate results
        src/filtering/s002_join_sbs_to_mbs_with_filtering.sh $DONT_KEEP_NESTED_SBS \
            -i $OUT_FILE.temp.isec.vcf.gz -o $OUT_FILE
    else
        trap "rm -rf $OUT_FILE.temp.vcf.gz $OUT_FILE.temp.vcf.gz.tbi" RETURN

        #Apply blacklist filter
        bcftools isec -c none -n=$(expr 1 + $INVERT_BLACKLIST) -w1 $OUT_FILE.temp.vcf.gz $BLACKLIST \
            -Oz -o $OUT_FILE

        tabix -p vcf $OUT_FILE
    fi
}

#Function for processing a patient VCF to CSV
function process_vcf() {
    #Input and output
    local IN_FILE=$1
    local OUT_FILE=$2

    #Filter variants with a blacklist if available
    if [[ $BLACKLIST != "" ]]; then
        filter_vcf_with_blacklist $IN_FILE "$OUT_FILE.blacklist_filtered.vcf.gz"
        IN_FILE="$OUT_FILE.blacklist_filtered.vcf.gz"
        trap "rm -rf $IN_FILE $IN_FILE.tbi" RETURN
    fi

    #Header
    echo -ne "CHROM\tPOS\tREF\tALT" > $OUT_FILE
    echo -ne $INFO_QUERY_STRING | sed 's/%//g' >> $OUT_FILE
    perl -e '
        @samples = split /\t/, $ARGV[0];
        @formats = split /\\t%/, $ARGV[1];
        shift @formats;

        foreach $sample (@samples) {
            foreach $format (@formats) {
                print "\t$sample.$format";
            }
        }

        print "\n";
    ' "$(bcftools view -h $IN_FILE | grep -v "^##" | cut -f10-)" $FORMAT_QUERY_STRING \
        >> $OUT_FILE

    #Data
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT$INFO_QUERY_STRING[$FORMAT_QUERY_STRING]\n" $IN_FILE \
        >> $OUT_FILE
}

#Function for simply copying the input to target
function copy_csv() {
    #Input and output
    local IN_FILE=$1
    local OUT_FILE=$2

    cp $IN_FILE $OUT_FILE
}

#Process all patient files
for FILE in $(ls $IN_DIR/*$IN_SUF* | sed 's/@$//' | grep -vP "\.tbi$"); do
    #Apply file name regex
    if [[ $(basename $FILE | grep $FILENAME_REGEX"" | wc -l) -gt 0 ]]; then
        #Simply copy if input is CSV
        if [[ $CSV == 0 ]]; then
            process_vcf $FILE $(new_file_name $FILE $IN_SUF $OUT_SUF) &
        else
            copy_csv $FILE $(new_file_name $FILE $IN_SUF $OUT_SUF) &
        fi
    fi
done

wait

