#!/bin/bash

#Extracts CHROM, POS and cluster number from subclonal cluster variant tables.

#Options
    #Defaults
    IN_DIR="data/raw/subclones/"
    IN_SUF="_cluster_variants.csv"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i | --input-directory )        shift
                                            IN_DIR=$1 ;;
            -s-suffix | --subclone-suffix ) shift
                                            IN_SUF=$1
        esac
        shift
    done

#Name variables
    OUT_DIR="data/features/$(basename $IN_DIR)/cluster_variants/"
    OUT_SUF=".tsv"

#Make sure target folder exist
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Function for finding column index given column name
get_column_index() {
    #Parameters
    local FILE=$1
    local COL_NAME=$2

    #Compute and return column index in STDOUT
    head -1 $FILE |
        perl -e '
            while(<STDIN>) {
                $fs_col = -1;
                chomp;
                @fields = split "\t";
                for $i (0 .. $#fields) {
                    $fs_col = $i + 1 if ($fields[$i] eq $ARGV[0]);
                }
            }
            print $fs_col;
        ' $COL_NAME
}

#Create the variant cluster tracks for each patient
for FILE in $IN_DIR/*$IN_SUF; do
    #Name variables
    NEW_FILE=${OUT_DIR}`perl -e 'print $1 if ($ARGV[0] =~ m/ \/ ([^\/]+) $ARGV[1] $ /x)' $FILE $IN_SUF`$OUT_SUF

    #Find mutation column and cluster column indices
    IND_MUTATION=$(get_column_index $FILE "mutation_id")
    IND_CLUSTER=$(get_column_index $FILE "cluster_id")

    #Output header
    echo -e "CHROM\tPOS\tcluster" > $NEW_FILE

    #Output data
    paste \
        <(cut -f$IND_MUTATION $FILE | tail -n +2 | perl -pe 's/^([^:]+):([^:]+).*\n$/$1\t$2\n/') \
        <(cut -f$IND_CLUSTER $FILE | tail -n +2) |
        perl -pe 's/^chrX/chr23/; s/^chrY/chr24/; s/^chrM/chr25/' |
        sort -k1,1V -k2,2n |
        perl -pe 's/^chr23/chrX/; s/^chr24/chrY/; s/^chr25/chrM/' |
        uniq \
            >> $NEW_FILE
done

