#!/bin/bash

#Script to get each SBS variant's transcriptional strand in terms of pyrimidine
#mutation

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
    IN_DIR="data/interim/filt_set_$NUM_ID/variants_sbs_ctx/"
    IN_SUF=".tsv"

    OUT_DIR="data/features/filt_set_$NUM_ID/variants_sbs/"
    OUT_SUF="_sbs_cts.tsv"

    GENES_TRANSCRIPTS="data/interim/genome_annotations/GRCh38.GENCODE.v29.known_genes_whole.bed"
    GENES_EXONS="data/interim/genome_annotations/GRCh38.GENCODE.v29.known_genes_exons.bed"
    CHR_FILE="data/interim/reference_genomes/GRCh38.chr_order.txt"

#Make sure target directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*$OUT_SUF

#Create #Patient-wise function to get transcriptional strand of mutation
function patient_sbs_ts() {
    #Variables
    local FILE=$1
    local IN_SUF=$2
    local OUT_SUF=$3
    local PATIENT=$(echo $FILE | sed s/$IN_SUF$// | sed 's,'^"$IN_DIR"',,')
    local NEW_FILE=${OUT_DIR}${PATIENT}$OUT_SUF

    #Header
    head -1 $FILE |
        awk '{$8 = $8 "\ttranscr_str\ttranscr_info"; print}' OFS='\t' \
        > $NEW_FILE

    #Primary categorisation of variants to intergenic [.] / transcribed single strand
    #[+/-] / transcr. both strands [=]
    tail -n +2 $FILE |
        #BED conversion
        awk '{$2 = ($2 - 1) "\t" $2; print}' OFS='\t' |
        #Intersect variants with whole genes
        bedtools intersect -c -sorted -g $CHR_FILE -a stdin -b $GENES_TRANSCRIPTS \
        > ${IN_DIR}${PATIENT}$IN_SUF.gene_intersect

    #Get variants not aligned (0) with transcripts and remove intersect count column
    perl -ane '
        if ($F[$#F] == 0) {
            $F[8] = "$F[8]\t.\t.";
            print(join("\t", @F[0 .. $#F - 1]) . "\n")
        }' ${IN_DIR}${PATIENT}$IN_SUF.gene_intersect \
        > ${IN_DIR}${PATIENT}$IN_SUF.gene_strand

    #Get variants aligned on both strands (2) and remove intersect count column
    perl -ane '
        if ($F[$#F] == 2) {
            $F[8] = "$F[8]\tboth\t.";
            print(join("\t", @F[0 .. $#F - 1]) . "\n")
        }' ${IN_DIR}${PATIENT}$IN_SUF.gene_intersect \
        >> ${IN_DIR}${PATIENT}$IN_SUF.gene_strand

    #Reintersect singly aligned variants (1) to get strand and exon/intron and
    #assign transcribed/untranscribed strand
    awk '$NF == 1' ${IN_DIR}${PATIENT}$IN_SUF.gene_intersect |
        #Intersect with whole genes for strand
        bedtools intersect -wa -wb -sorted -g $CHR_FILE -a stdin -b $GENES_TRANSCRIPTS |
        #Intersect with exons for exon / -1 (aka intron)
        bedtools intersect -loj -sorted -g $CHR_FILE -a stdin -b $GENES_EXONS |
        #Replace -1 with intron and reposition gene strand and exon/intron columns
        perl -ane '
            if ($F[$#F - 1] == -1) {$F[$#F - 1] = "intron"}
            print (join("\t", (@F[0 .. 8], $F[$#F - 6], $F[$#F - 1], @F[9 .. $#F - 13])) . "\n")
        ' |
        #Compute transcriptional strand for each variant
        awk '{if($9 == $10) {$10 = "transcribed"} else {$10 = "untranscribed"} print}' OFS='\t' \
        >> ${IN_DIR}${PATIENT}$IN_SUF.gene_strand

    #Sort and convert back to tsv (single position) output final to target file
    sed 's/^chrM/chrZ/' ${IN_DIR}${PATIENT}$IN_SUF.gene_strand |
        sort -k1,1V -k2,2n |
        sed 's/^chrZ/chrM/' |
        cut -f1,3- \
        >> $NEW_FILE

    #Remove temporary file
    rm ${IN_DIR}${PATIENT}$IN_SUF.gene_intersect
    rm ${IN_DIR}${PATIENT}$IN_SUF.gene_strand
}

#Loop over patients
for FILE in $IN_DIR/*$IN_SUF; do
    #SBS with DBS excluded
    patient_sbs_ts $FILE $IN_SUF $OUT_SUF &
done

wait

