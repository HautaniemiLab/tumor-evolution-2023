#!/bin/bash

#Mini shell workflow for preparing resources for signature analysis. This
#includes:
# - Acquiring genome fastas for different genome builds
# - Computing dinucleotide and trinucleotide frequencies
# - Downloading GENCODE
# - Extracting known exons/transcripts from GENCODE
# - Downloading reference signatures (only COSMIC is automatic)
# - Correcting signatures to GRCh38 trinucleotide frequencies and as TSV

#Has an option to use alternative main genome (user supplied)
#If targets are supplied, reference signatures for them are also computed
#using assumed or given genome

#TODO: default genome in some config instead?

#Options
    #Defaults
    NOT_DEFAULT_GENOME=""
    NO_MAIN_GENOME=""
    MAIN_GENOME="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    MAIN_CHR_FILTER="M|Y|_|EBV|^[HKS]"
    TARGETS=""

    #Command line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -g | --genome-fasta )   shift
                                    MAIN_GENOME=$1
                                    NOT_DEFAULT_GENOME="TRUE" ;;
            --chromosome-filter )   shift
                                    MAIN_CHR_FILTER=$1 ;;
            --no-main-genome )      NO_MAIN_GENOME="TRUE" ;;
            -t | --targets )        shift
                                    TARGETS=$1
        esac
        shift
    done

#Get genome fastas
    #GRCh38d1vd1 (default)
    if [[ $NO_MAIN_GENOME == "" && $NOT_DEFAULT_GENOME == "" ]]; then
        echo "Getting GRCh38d1vd1 genome fasta"
        src/data/s000_get_GRCh38d1vd1_fasta_link.sh
    fi

    #hs37d5
    echo "Getting hs37d5 genome fasta"
    src/data/s010_download_hs37d5_fasta.sh

    #hg19
    echo "Getting hg19 genome fasta"
    src/data/s020_download_hg19_fasta.sh

#Compute genome and target i.e. exome 2- and 3-mer counts
    echo "Computing genome 2-mer and 3-mer counts"
    if [[ $NO_MAIN_GENOME == "" ]]; then
        src/features/s100_genome_dinuc_trinuc_counts.sh -g "$MAIN_GENOME" --chromosome-filter "$MAIN_CHR_FILTER"
    fi
    src/features/s100_genome_dinuc_trinuc_counts.sh -g "data/raw/reference_genomes/hs37d5.fa" --chromosome-filter "^(M|GL|NC|hs37d5)"
    src/features/s100_genome_dinuc_trinuc_counts.sh -g "data/raw/reference_genomes/hg19.fa" --chromosome-filter "M|_"

#Download and prepare GENCODE for GRCh38
    #Download raw GENCODE
    echo "Getting GENCODE"
    src/data/s100_download_gencode.sh

    #Get known transcripts and exons from GENCODE
    echo "Extracting GENCODE known transcripts and exons"
    src/features/s000_sbs_prep_gene_annotations.sh

#Download and prepare signatures
    #Download COSMIC v2 signatures
    echo "Getting COSMIC v2 signatures"
    src/data/s200_download_cosmic_signatures.sh

    #Download COSMIC v3 signatures and convert them to TSV
    echo "Getting COSMIC v3 signatures"
    src/data/s210_prepare_sigProfiler_signatures.sh

    #Correct signatures based on genome kmer counts
    if [[ $NO_MAIN_GENOME == "" ]]; then
        echo "Correcting signatures to GRCh38 kmer counts"
        src/features/s110_correct_ref_sigs.sh -t "$MAIN_GENOME" --relabel-indels
    fi

#Handle targets
    if [[ $TARGETS != "" ]]; then
        src/s010_prepare_signatures_for_custom_genome.sh -g "$MAIN_GENOME" -t "$TARGETS"
    fi

