#!/bin/bash

#Script to produce two tables for matching original nearby (gap < 10bp)
#phased mutations with their merged counterparts for dNdScv

REFERENCE="/path/to/GRCh38.d1.vd1.fa"

#Work directory
OUT_DIR="phased_mutations"
mkdir -p $OUT_DIR

#Function to process one patient's variants to merge and gaps between
function process_phased_mutations_and_gaps() {
    local PATIENT=$1
    local VCF=$2
    local PREFIX="$3"

    #Look only at variants with phasing information (PGT:PID)
    bcftools view -H $VCF |
        cut -f-2,4,5,9,10 |
        grep -P ":PGT:PID:" |
        perl -e '
            use warnings;
            use strict;

            sub which_index {
                my ($col, $field_name) = @_;

                my @col_split = split /:/, $col;

                for my $i (0 .. $#col_split) {
                    return ($i) if ($col_split[$i] eq $field_name);
                }

                return -1;
            }

            sub get_index {
                my ($col, $index) = @_;

                return (split /:/, $col)[$index];
            }

            print "chr\tpos1\tpos2\tref1\tref2\tmut1\tmut2\tindex\tgap\tgapChr\tgapStart\tgapEnd\n";

            my $index = 1;
            my ($prev_chr, $prev_pos, $prev_ref, $prev_mut, $prev_pgt, $prev_pid) = ("", "-inf", "", "", "", "");
            my @prev_line;

            while (<STDIN>) {
                chomp;
                my @line = split /\t/;

                my ($pgt_ind, $pid_ind) = (which_index($line[4], "PGT"), which_index($line[4], "PID"));
                my ($chr, $pos, $ref, $mut, $pgt, $pid) = (@line[0 .. 3], get_index($line[5], $pgt_ind), get_index($line[5], $pid_ind));

                #Same chromosome, PGT, PID and distance is within 10
                my $gap = $pos - $prev_pos - length($prev_ref);

                #gap < 10 to better consider indels; gap >= 0 to avoid first deletion covering second variant
                if (
                    $prev_chr eq $chr and
                    $gap < 10 and
                    $gap >= 0 and
                    $prev_pgt eq $pgt and
                    $prev_pid eq $pid
                ) {
                    print join("\t", ($prev_chr, $prev_pos, $pos, $prev_ref, $ref, $prev_mut, $mut, $index++)) .
                        "\t$gap" . "\t$chr\t" . ($pos - 1 - $gap) . "\t" . ($pos - 1) . "\n";
                }

                ($prev_chr, $prev_pos, $prev_ref, $prev_mut, $prev_pgt, $prev_pid) = ($chr, $pos, $ref, $mut, $pgt, $pid);
                @prev_line = @line;
            }
        ' \
            > $OUT_DIR/${PREFIX}vars_to_merge/$PATIENT.tsv

    awk '{if (NR > 1 && $9 > 0) {print $10, $11, $12}}' OFS="\t" $OUT_DIR/${PREFIX}vars_to_merge/$PATIENT.tsv |
        bedtools getfasta -fi $REFERENCE -bed stdin -s -tab -fo stdout |
        sed 's/()//; s/[:\-]/\t/g' |
        cat <(echo -e "gapChr\tgapStart\tgapEnd\tgapSeq") - \
            > $OUT_DIR/${PREFIX}vars_to_merge.gaps/$PATIENT.tsv
}

#Main processing function
function generate_phased_mutation_tables() {
    local VCF_FILES_TBL=$1
    local PREFIX="$2"

    mkdir -p $OUT_DIR/${PREFIX}vars_to_merge $OUT_DIR/${PREFIX}vars_to_merge.gaps

    #Determine variants (always at most two) to merge
    while read PATIENT VCF; do
        process_phased_mutations_and_gaps $PATIENT $VCF "$PREFIX" &
    done < <(cut -f1,3 $VCF_FILES_TBL | tail -n +2)

    wait

    #Merge results
    Rscript -e '
        suppressMessages(library(tidyverse))

        options(readr.show_col_types=FALSE)

        args = commandArgs(trailingOnly=TRUE)

        in_tbl = args[1]
        prefix = paste0("phased_mutations/", args[2])

        samples = read_tsv(in_tbl)

        mut_tables = lapply(
            samples$patient,
            function(patient) {
                vars_to_merge = read_tsv(paste0(prefix, "vars_to_merge/", patient, ".tsv"))
                vars_to_merge_gaps = read_tsv(paste0(prefix, "vars_to_merge.gaps/", patient, ".tsv"))

                vars_processed = vars_to_merge %>%
                    merge(vars_to_merge_gaps, by = c("gapChr", "gapStart", "gapEnd"), all.x=T) %>%
                    mutate(gapSeq = ifelse(gap == 0, "", gapSeq))

                list(
                    "split" = cbind(
                        "patient" = rep(patient, nrow(vars_processed)),
                        rbind(
                            vars_processed %>% select(chr, pos=pos1, ref=ref1, mut=mut1, index),
                            vars_processed %>% select(chr, pos=pos2, ref=ref2, mut=mut2, index)
                        ) %>%
                            arrange(chr, pos)
                    ),
                    "merged" = vars_processed %>%
                        mutate(
                            patient = rep(patient, nrow(vars_processed)),
                            ref = paste0(ref1, gapSeq, ref2),
                            mut = paste0(mut1, gapSeq, mut2)
                        ) %>%
                        select(patient, chr, pos=pos1, mut, ref, index)
                )
            }
        )

        lapply(mut_tables, function(x) x$split) %>%
            do.call(rbind, .) %>%
            arrange(patient, index) %>%
            write_tsv(paste0(prefix, "phased_mutations.split.tsv"))
        lapply(mut_tables, function(x) x$merged) %>%
            do.call(rbind, .) %>%
            arrange(patient, index) %>%
            write_tsv(paste0(prefix, "phased_mutations.merged.tsv"))
    ' $VCF_FILES_TBL "$PREFIX"
}

#R scripts to produce a table of VCF files to process
#Discovery
Rscript -e '
    suppressMessages(library(tidyverse))

    suppressMessages(read_tsv("data/samples.tsv")) %>%
        filter(gistic) %>%
        select(patient, sample) %>%
        mutate(vcf = paste0("data/variants_vcf/", patient, ".vcf.gz")) %>%
        write_tsv("phased_mutations/vcf_files.tsv")
'

#Validation
Rscript -e '
    suppressMessages(library(tidyverse))

    suppressMessages(read_tsv("data/samples.validation.tsv")) %>%
        #Purity >= 10% non-duplicate tumour
        filter(purity >= 0.1, isMain, !normal) %>%
        #Choose sample with the same GISTIC prioritisation
        mutate(
            psNonAsc = !(tissueType %in% c("ascites", "pleural fluid")),
            psNonAdn = !(tissueType %in% c("adnex", "fallopian tube", "ovary")),
            psNaive = sampleTime == "primary",
            priorityScore = 4 * psNonAsc + 2 * psNonAdn + 1 * psNaive
        ) %>%
        group_by(patient) %>%
        arrange(priorityScore, purity) %>%
        slice_tail(n = 1) %>%
        ungroup() %>%
        #Get VCF path
        select(patient, sample) %>%
        mutate(vcf = paste0("data/variants_vcf.validation/", patient, ".vcf.gz")) %>%
        write_tsv("phased_mutations/validation.vcf_files.tsv")
'

#Generate the tables
generate_phased_mutation_tables "phased_mutations/vcf_files.tsv" ""
generate_phased_mutation_tables "phased_mutations/validation.vcf_files.tsv" "validation."
