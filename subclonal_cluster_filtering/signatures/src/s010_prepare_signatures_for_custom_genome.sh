#!/bin/bash

#Mini shell workflow for adjusting signatures for custom genome builds or
#regions. In other words
# - Computes dinucleotide and trinucleotide frequencies
# - Corrects signatures to GRCh38 trinucleotide frequencies and as TSV

#Options
    #Defaults
    GENOME="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    TARGETS=""
    CHR_FILTER="M|Y|_|EBV|^[HKS]"

    while [[ $# -gt 0 ]]; do
        case $1 in
            -g | --genome )         shift
                                    GENOME=$1 ;;
            -t | --targets )        shift
                                    TARGETS=$1 ;;
            --chromosome-filter )   shift
                                    CHR_FILTER=$1
        esac
        shift
    done

#Compute genome and/or target (e.g. exome) 2- and 3-mer counts
    if [[ $GENOME != "data/raw/reference_genomes/GRCh38.d1.vd1.fa" ]]; then
        echo "Computing $GENOME 2-mer and 3-mer counts"
        src/features/s100_genome_dinuc_trinuc_counts.sh -g $GENOME --chromosome-filter "$CHR_FILTER"
    fi

    if [[ $TARGETS != "" ]]; then
        echo "Computing $TARGETS 2-mer and 3-mer counts"
        src/features/s100_targets_dinuc_trinuc_counts.sh -g $GENOME -t $TARGETS --chromosome-filter "$CHR_FILTER"
    fi

#Adjust signatures
    #Correct signatures based on genome kmer counts
    if [[ $GENOME != "data/raw/reference_genomes/GRCh38.d1.vd1.fa" ]]; then
        echo "Correcting signatures to $GENOME kmer counts"
        src/features/s110_correct_ref_sigs.sh -t $GENOME
    fi

    #Correct signatures based on target kmer counts
    if [[ $TARGETS != "" ]]; then
        echo "Correcting signatures to $TARGETS kmer counts"
        src/features/s110_correct_ref_sigs.sh -t $TARGETS
    fi

