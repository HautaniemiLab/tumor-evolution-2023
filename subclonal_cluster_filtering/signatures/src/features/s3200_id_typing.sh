#!/bin/bash

#Script to type InDels i.e. assign mutation type and to same classification as in
#Alexandrov et al. ? (bioRxiv 2018)

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
    TYPING_SCRIPT="src/features/s3201_id_typing.pl"

    IN_DIR="data/interim/filt_set_$NUM_ID/variants_id/"
    IN_SUF=".tsv"
    
    OUT_DIR="data/features/filt_set_$NUM_ID/variants_id/"
    OUT_SUF=".tsv"

#Mke sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Patient-wise function to get mutation type
function patient_id_types() {
    #Variables
    local FILE=$1
    local PATIENT=$(echo $FILE | sed s/$IN_SUF$// | sed 's,'^"$IN_DIR"',,')
    local NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Header
    head -1 $FILE |
        awk '{$5 = $5 "\ttype"; print}' OFS='\t' \
        > $NEW_FILE

    #Get In/Del and length for mutation type categorisation and output as a temp bed
    #file
    tail -n +2 $FILE |
        perl -ane '
            if (length $F[2] == 1) {$F[4] = "$F[4]\tINS"} else {$F[4] = "$F[4]\tDEL"}
            $F[4] = "$F[4]\t" . abs(length($F[3]) - length($F[2]));
            $F[1] = ($F[1] - 1) . "\t$F[1]";
            print(join("\t", @F) . "\n");
        ' \
        > ${IN_DIR}${PATIENT}_id_temp.bed

    #Get flanking sequence for typing:
    #Get positions and length
    cut -f-3,8 ${IN_DIR}${PATIENT}_id_temp.bed |
        #Expand positions, keep original 1-based positions in the strand column
        awk '{len = $4; $4 = ".\t.\t" $3; $2 = $2 - len + 2; $3 = $3 + 6 * len; print}' OFS='\t' |
        #Extract respective sequences from reference genome
        bedtools getfasta -fi $REF_FASTA -bed stdin -s -tab -fo stdout |
        #Decompose getfasta output to tab-separated format, extract 'strand' field as
        #position
        perl -pe 's/ (.+) : ([0-9]+) - ([0-9]+) \( ([0-9]+) \) (.+) /$1\t$4$5/x' |
        #Convert single point position back to BED
        awk '{$2 = ($2 - 1) "\t" $2; print}' OFS='\t' |
        #Intersect sequences with ID variants
        bedtools intersect -a stdin -b ${IN_DIR}${PATIENT}_id_temp.bed -wa -wb |
        #Remove duplicate position columns
        cut -f4,5,7- |
        #Reposition context column
        perl -ane 'print (join("\t", (@F[1 .. 7], $F[0], @F[8 .. $#F])) . "\n")' |
        #Run indel typing script
        $TYPING_SCRIPT |
        #Remove excess columns used for typing
        cut -f-6,10- \
        >> $NEW_FILE

    #Clear temp file
    rm ${IN_DIR}${PATIENT}_id_temp.bed
}

#Loop over files to determine strand and primary mutation type, then extract from
#reference fasta its context
for FILE in $IN_DIR/*$IN_SUF; do
    patient_id_types $FILE &
done

wait

