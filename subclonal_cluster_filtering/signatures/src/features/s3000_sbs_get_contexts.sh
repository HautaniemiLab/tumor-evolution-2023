#!/bin/bash

#Script to type SBS variants i.e. assign pyrimidine mutation type and strand as
#well as its trinucleotide context.

#Options
    #Defaults
    NUM_ID="0"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )     shift
                                NUM_ID=$1
        esac
        shift
    done

#Name variables
    REF_FASTA="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    IN_DIR="data/interim/filt_set_$NUM_ID/variants_sbs/"
    IN_SUF=".tsv"

    OUT_DIR="data/interim/filt_set_$NUM_ID/variants_sbs_ctx/"
    OUT_SUF=".tsv"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Patient-wise function to get mutation type, context and strand
function patient_sbs_contexts() {
    #Input variables
    local FILE=$1
    local IN_DIR=$2
    local IN_SUF=$3
    local OUT_DIR=$4
    local OUT_SUF=$5

    #Name variables
    local PATIENT=$(echo $FILE | sed s/$IN_SUF$// | sed 's,'^"$IN_DIR"',,')
    local NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Header
    head -1 $FILE |
        awk '{$5 = $5 "\ttype\tsubtype\tstrand"; print}' OFS='\t' \
        > $NEW_FILE

    #Get type and strand as new columns, and output as a temp bed file
    tail -n +2 $FILE |
        perl -ane '
            $type = "$F[2]>$F[3]";
            if ($F[2] =~ /[CT]/) {
                $F[4] = "$F[4]\t$type\t+";
            } else {
                $type =~ tr/ACGT/TGCA/;
                $F[4] = "$F[4]\t$type\t-";
            }
            $F[1] = ($F[1] - 1) . "\t" . $F[1];
            print(join("\t", @F) . "\n");
        ' \
        > ${IN_DIR}${PATIENT}$IN_SUF.temp.bed

    #Finish up typing by getting contexts:
    #Convert temp file to bed with strand info and trinucleotide context positions,
    #then use bedtools getfasta to get contexts. Finally convert back to single SBS
    #positions and intersect with the temp file, reorder columns as final output
    #Get positions, type and strand
    cut -f-3,7,8 ${IN_DIR}${PATIENT}$IN_SUF.temp.bed |
        #Expand positions to trinucleotide, make strand sixth column
        awk '{$2--; $3++; $4 = $4 "\t1"; print}' OFS='\t' |
        #Extract trinucleotides from reference sequence
        bedtools getfasta -fi $REF_FASTA -bed stdin -s -tab -fo stdout |
        #Decompose getfasta output to tab-separated format
        perl -pe 's/ (.+) : ([0-9]+) - ([0-9]+) \( ([+-]) \) (.+) /$1\t$2\t$3\t$4$5/x' |
        #Convert back to single point positions
        awk '{$2++; $3--; print}' OFS='\t' |
        #Intersect trinucleotides with SBS variants
        bedtools intersect -a stdin -b ${IN_DIR}${PATIENT}$IN_SUF.temp.bed -wa -wb |
        #Remove duplicate position columns
        cut -f5,6,8- |
        #Reposition context column
        perl -ane 'print (join("\t", (@F[1 .. 6], $F[0], @F[7 .. $#F])) . "\n")' \
        >> $NEW_FILE

    #Clear temp file
    rm ${IN_DIR}${PATIENT}$IN_SUF.temp.bed
}

#Loop over files to determine strand and primary mutation type, then extract its
#context from reference fasta
for FILE in $IN_DIR/*$IN_SUF; do
    patient_sbs_contexts $FILE $IN_DIR $IN_SUF $OUT_DIR $OUT_SUF &
done

wait

