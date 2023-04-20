#!/bin/bash

#Counts the di- and trinucleotide frequencies in specified genome build (fasta)
#with Alexandrov et al. 2020 reference dinucleotides and central pyrimidines in
#trinucleotides. Alexandrov et al. 2020 (sigProfiler aka COSMIC v3) uses hs37d5
#while for COSMIC v2 GRCh37/hg19 should be used.
#Uses customisable grep filter to remove undesired chromosomes such as decoys

#Options
    #Defaults
    IN_FILE="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    CHR_FILTER="M|Y|_|EBV|^[HKS]"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -g | --genome-fasta )   shift
                                    IN_FILE=$1 ;;
            --chromosome-filter )   shift
                                    CHR_FILTER=$1
        esac
        shift
    done

#Name variables
    CHR_EXTRACTION_SCRIPT="src/features/s101_extract_chrs.pl"
    KMER_COUNTING_SCRIPT="src/features/s102_genome_kmer_count.pl"
    TBL_COMBINING_SCRIPT="src/features/s103_combine_kmer_count_tables.pl"

    GENOME=`echo $IN_FILE | perl -pe 's/.*\///; s/\..*//'`

    OUT_DIR="data/interim/reference_genomes/"
    CHR_DIR="${OUT_DIR}${GENOME}.chromosomes/"
    CHR_LIST="${CHR_DIR}chr_list"
    CHR_COUNTS_LIST="${CHR_DIR}chr_counts_list"

    DINUC_OUT_SUF=".dinuc_freq.tsv"
    TRINUC_OUT_SUF=".trinuc_freq.tsv"

#Make sure chromosome directory exists
    mkdir -p $CHR_DIR

#Create the list of chromosomes
    grep -P "^>" $IN_FILE | tr -d ">" | cut -f1 -d' ' \
        > $CHR_LIST

#Create individual one-string chromosome files
    $CHR_EXTRACTION_SCRIPT -l $CHR_LIST -b $IN_FILE -o $CHR_DIR

#Count triplets and doublets by chromosome with filter applied
    rm -f $CHR_COUNTS_LIST.*

    while read CHR; do
        #Triplets
        $KMER_COUNTING_SCRIPT -c ${CHR_DIR}$CHR -p 2 \
            > ${CHR_DIR}${CHR}.3-mer_counts.tsv &
        echo "${CHR_DIR}${CHR}.3-mer_counts.tsv" \
            >> ${CHR_COUNTS_LIST}.3
    done < <(grep -vP $CHR_FILTER $CHR_LIST)

wait

    while read CHR; do
        #Doublets
        $KMER_COUNTING_SCRIPT -c ${CHR_DIR}$CHR -k 2 -f \
            > ${CHR_DIR}${CHR}.2-mer_counts.tsv &
        echo "${CHR_DIR}${CHR}.2-mer_counts.tsv" \
            >> ${CHR_COUNTS_LIST}.2
    done < <(grep -vP $CHR_FILTER $CHR_LIST)

wait

#Combine tables
$TBL_COMBINING_SCRIPT -l ${CHR_COUNTS_LIST}.2 \
    > ${OUT_DIR}${GENOME}$DINUC_OUT_SUF
$TBL_COMBINING_SCRIPT -l ${CHR_COUNTS_LIST}.3 \
    > ${OUT_DIR}${GENOME}$TRINUC_OUT_SUF

