#!/bin/bash

#Shell script for splitting multiple base substitions (MBS i.e. MNPs) in biallelic
#VCF to SBS (i.e. SNPs) while keeping the original records. This is done to help
#filter DBS with high germline risk for at least one of the nested SBS. Adds a
#unique ID to the DBS and SBS to link them for later rejoining and filtering by
#using a custom INFO field MBS_ID. Outputs as VCF.GZ.

#Parameters
    #Defaults
    INPUT_FILE=$1
    OUTPUT_FILE=$2
    TEMP_DIR="$(pwd)/tmp.$$/"

    #Command Line
    while [ "$1" != "" ]; do
        case $1 in
            -i | --input )                  shift
                                            INPUT_FILE=$1 ;;
            -o | --output )                 shift
                                            OUTPUT_FILE=$1 ;;
            -t | --temporary-directory )    shift
                                            TEMP_DIR=$1
        esac
        shift
    done

#Make sure inputs and output are given
if [[ $INPUT_FILE == "" ]]; then
    echo "Input VCF must be given as second argument or via -i or --input" >&2
    exit 1;
fi
if [[ $OUTPUT_FILE == "" ]]; then
    echo "Output File must be given as third argument or via -o or --output" >&2
    exit 1;
fi

#Make sure input files exist
if [[ ! -f $INPUT_FILE ]]; then
    echo "Input VCF $INPUT_FILE doesn't exist!" >&2
    exit 1;
fi

#Create temporary directory if it doesn't exist
if [[ ! -d $TEMP_DIR ]]; then
    mkdir $TEMP_DIR
    trap "rm -rf $TEMP_DIR" EXIT
fi

#Print all SBS within MBS add running number identifiers to QUALs. Output as VCF.GZ
bcftools view -H $INPUT_FILE |
    perl -ane '
        BEGIN{$num = 1}
        if (length($F[3]) == length($F[4]) and length($F[3]) > 1) {
            print join("\t", @F[0 .. 6]) . "\tMBS_ID=$num;" . join("\t", @F[7 .. $#F]) . "\n";

            for ($i=0; $i<length($F[3]); $i++) {
                $pos = $F[1] + $i;
                $ref = substr($F[3], $i, 1);
                $alt = substr($F[4], $i, 1);
                $id = "$num." . ($i < 9 ? "0" : "") . ($i + 1);

                print "$F[0]\t$pos\t$F[2]\t$ref\t$alt\t" . join("\t", @F[5 .. 6]) . "\tMBS_ID=$id;" . join("\t", @F[7 .. $#F]) . "\n";
            }

            $num++
        } else {
            print;
        }
    ' |
    #Sort the records by coordinates
    perl -pe 's/^chrX/chr23/; s/^chrY/chr24/; s/^chrM/chr25/' |
    sort -k1,1V -k2,2V -k4,4V -k5,5V -T $TEMP_DIR |
    perl -pe 's/^chr23/chrX/; s/^chr24/chrY/; s/^chr25/chrM/' |
    #Add header
    cat <(bcftools view -h $INPUT_FILE | perl -ne 'if (not m/^##/) {print "##INFO=<ID=MBS_ID,Number=1,Type=String,Description=\"MBS identifier\">\n"}; print') - |
    #Output as VCF.GZ
    bcftools view -Oz -o $OUTPUT_FILE

#Index the output VCF.GZ
tabix -p vcf $OUTPUT_FILE

