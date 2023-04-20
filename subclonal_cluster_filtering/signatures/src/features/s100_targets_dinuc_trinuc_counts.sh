#!/bin/bash

#Counts the di- and trinucleotide frequencies in specified targeted regions of
#genome build (fasta), e.g. exomes, with Alexandrov et al. 2020 reference
#dinucleotides and central pyrimidines in trinucleotides. Alexandrov et al. 2020
#(sigProfiler aka COSMIC v3)uses hs37d5 while for COSMIC v2 GRCh37/hg19 should
#be used.
#Uses customisable grep filter to remove undesired chromosomes such as decoys

#Options
    #Defaults
    GENOME_FILE="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    OVERWRITE_GENOME=""
    CHR_FILTER="M|Y|_|EBV|^[HKS]"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -g | --genome-fasta )   shift
                                    GENOME_FILE=$1 ;;
            --overwrite-genome )    OVERWRITE_GENOME="TRUE" ;;
            --chromosome-filter )   shift
                                    CHR_FILTER=$1 ;;
            -t | --targets )        shift
                                    TARGETS_FILE=$1
        esac
        shift
    done

#Make sure targets is given and exists
if [[ ! -f $TARGETS_FILE ]]; then
    echo "Supply a valid targets (BED) file!" >&2
    exit 1;
fi

#Name variables
    CHR_EXTRACTION_SCRIPT="src/features/s101_extract_chrs.pl"
    KMER_COUNTING_SCRIPT="src/features/s102_genome_kmer_count.pl"

    GENOME=`echo $GENOME_FILE | perl -pe 's/.*\///; s/\..*//'`
    TARGET=`echo $TARGETS_FILE | perl -pe 's/.*\///; s/\..*//'`

    OUT_DIR="data/interim/reference_genomes/"
    CHR_DIR="${OUT_DIR}${GENOME}.chromosomes/"
    CHR_LIST="${CHR_DIR}chr_list"
    INTERVALS="${OUT_DIR}${TARGET}.bed"

    DINUC_OUT_SUF=".dinuc_freq.tsv"
    TRINUC_OUT_SUF=".trinuc_freq.tsv"

#Check if the genome files are already present or if it should be overwritten
    if [[ $OVERWRITE_GENOME != "" || ! -d $CHR_DIR || ! -f $CHR_LIST ]]; then
        #Make sure chromosome directory exists
        mkdir -p $CHR_DIR

        #Create the list of chromosomes:
        grep -P "^>" $GENOME_FILE | tr -d ">" | cut -f1 -d' ' \
            > $CHR_LIST

        #Create individual one-string chromosome files
        $CHR_EXTRACTION_SCRIPT -l $CHR_LIST -b $GENOME_FILE -o $CHR_DIR
    fi

#Prepare targets file
    bedtools merge -i $TARGETS_FILE |
        #apply chromosome filter
        perl -e '
            while(<STDIN>) {
                @F = split /\t/;
                print if ($F[0] !~ m/$ARGV[0]/);
            }
        ' "$CHR_FILTER" \
            > $INTERVALS

#Count all triplets and doublets in the intervals
    $KMER_COUNTING_SCRIPT -d ${CHR_DIR} -l $CHR_LIST -intervals $INTERVALS -p 2 -padding 1 \
        > ${OUT_DIR}${TARGET}$TRINUC_OUT_SUF &

    $KMER_COUNTING_SCRIPT -d ${CHR_DIR} -l $CHR_LIST -intervals $INTERVALS -k 2 -f \
        > ${OUT_DIR}${TARGET}$DINUC_OUT_SUF &

wait

