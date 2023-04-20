#!/bin/bash

#Shell script for reconstructing multiple base substitutions by its nested single
#base substitutions based on their existence as individual records.Potentially
#joins passing SBS to smaller MBS. Adds phasing FORMAT fields PGT, PID and PS to
#separated SBS due to filtering. Outputs as VCF.GZ.

#Filtering and printing rules:
# - print nothing if all SBS are filtered
# - print MBS if MBS and all of its SBS are unfiltered (SBS records discarded)
# - print all individual SBS if MBS unfiltered but at least one SBS filtered
#   - contiguous SBS are combined to MBS
#   - This can be disabled via an argument

#Parameters
    #Defaults
    INPUT_FILE=$1
    OUTPUT_FILE=$2
    DONT_KEEP_NESTED_SBS="0"
    TEMP_DIR="$(pwd)/tmp.$$/"

    #Command Line
    while [ "$1" != "" ]; do
        case $1 in
            -i | --input )                  shift
                                            INPUT_FILE=$1 ;;
            -o | --output )                 shift
                                            OUTPUT_FILE=$1 ;;
            --dont-keep-nested-sbs )        DONT_KEEP_NESTED_SBS="1" ;;
            -t | --temporary-directory )    shift
                                            TEMP_DIR=$1
        esac
        shift
    done

#Make sure inputs and output are given
if [[ $INPUT_FILE == "" ]]; then
    echo "Input VCF must be given as first argument or via -i or --input" >&2
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

#Process the MBS records to have either MBS or SBS representation
bcftools view -H $INPUT_FILE |
    perl -ane 'if ($F[7] =~ m/^MBS_ID=(\d+\.?\d{0,2});/) {print "$1"} print "\t$_"' |
    sort -V |
    cut -f2- |
    perl -e '
        use strict;
        use warnings;
        use List::Util qw(max);

        my $dont_keep_nested_sbs = $ARGV[0];

        my @F;
        my $mbs_id = -1;
        my $mbs_line;
        my @sbs_lines = ();
        my @sbs_ids = ();

        #Loop over records
        while (<STDIN>) {
            @F = split /\t/;

            #Distinguish MBS (new) and SBS (same MBS) lines
            #New MBS
            if ($F[7] =~ m/^MBS_ID=(\d+);/) {
                #Handle SBS or MBS printing if MBS was not already printed; skip for first line
                filter_and_print_record($mbs_line, \@sbs_lines, \@sbs_ids) if ($mbs_line);

                #Initialise new MBS
                $mbs_line = $_;
                $mbs_id = $1;
                @sbs_lines = ();
                @sbs_ids = ();
            } elsif ($F[7] =~ m/^MBS_ID=(\d+)\.(\d+);/) {
                #Add SBS belonging to MBS; orphan SBSs are discarded
                if ($1 == $mbs_id) {
                    push @sbs_lines, $_;
                    push @sbs_ids, $2;
                }
            } else {
                #Print non-MBS line
                print;
            }
        }
        filter_and_print_record($mbs_line, \@sbs_lines, \@sbs_ids) if ($mbs_line);

        #Subroutine for deciding whether to print MBS or SBS
        sub filter_and_print_record {
            my ($mbs_line, $sbs_lines_ref, $sbs_ids_ref) = @_;

            #Expected number of SBS
            my $expected_sbs_count = length((split(/\t/, $mbs_line))[3]);

            #Print MBS if all SBS have been retained
            if ($expected_sbs_count == scalar @$sbs_ids_ref) {
                clean_info_and_print_record($mbs_line);
                return;
            }

            #Print nothing if no SBS has been retained or their printing has been disabled
            if (scalar @$sbs_ids_ref == 0 or $dont_keep_nested_sbs) {
                return;
            }

            #Print individual SBS; if any neighbouring SBS are kept, join them as MBS
            my @prev_record;
            my $prev_id = -1;
            my @curr_record;
            my $curr_id;
            my ($pid, $ps);

            for my $i (0 .. $#{$sbs_ids_ref}) {
                #Read new SBS record
                @curr_record = split /\t/, $sbs_lines_ref->[$i];
                $curr_id = $sbs_ids_ref->[$i];

                #Is there a gap between SNPs
                if ($prev_id > 0 and $curr_id == $prev_id + 1) {
                    $prev_record[3] .= $curr_record[3];
                    $prev_record[4] .= $curr_record[4];
                } else {
                    #Handle previous record if it was not yet printed
                    ($pid, $ps) = phase_and_print_record(\@prev_record, $pid, $ps) if ($prev_id > 0);

                    @prev_record = @curr_record;
                }

                $prev_id = $curr_id;
            }

            #Print the trailing record
            phase_and_print_record(\@prev_record, $pid, $ps);
        }

        #Subroutine for cleaning MBS_ID INFO annotation and printing record
        sub clean_info_and_print_record {
            my $record = shift;

            my @fields = split /\t/, $record;
            $fields[7] =~ s/^MBS_ID=\d+\.?\d*;//;

            print join("\t", @fields);
        }

        #Subroutine for phasing and printing record
        sub phase_and_print_record {
            my ($record_ref, $pid, $ps) = @_;

            #Make sure phasing fields are defined
            ($pid, $ps) = define_pid_ps($record_ref) if (not $pid);

            #Update the phasing fields
            update_formats($record_ref, $pid, $ps);

            #Print record
            clean_info_and_print_record(join("\t", @{$record_ref}));

            return ($pid, $ps);
        }

        #Subroutine for defining FORMAT fields PID and PS
        sub define_pid_ps {
            my $record_ref = shift;

            #Define PID, and PS using POS, REF, ALT
            return ("$record_ref->[1]_$record_ref->[3]_$record_ref->[4]", $record_ref->[1]);
        }

        #Subroutine for updating FORMAT fields GT, PGT, PID, and PS if missing phasing information
        sub update_formats {
            my ($record_ref, $pid, $ps) = @_;

            my @formats = split /:/, $record_ref->[8];
            my ($gt_i, $pgt_i, $pid_i, $ps_i) = (-1, -1, -1, -1);

            #Look for indices for the fields
            for my $i (0 .. $#formats) {
                $gt_i = $i if($formats[$i] eq "GT");
                $pgt_i = $i if($formats[$i] eq "PGT");
                $pid_i = $i if($formats[$i] eq "PID");
                $ps_i = $i if($formats[$i] eq "PS");
            }

            #Insert PGT, PID, and PS if they are all missing
            if ($pgt_i < 0 and $pid_i < 0 and $ps_i < 0) {
                $formats[5] .= ":PGT:PID:PS";
                $record_ref->[8] = join(":", @formats);
                $record_ref->[8] =~ s/:+/:/g;
                my $inserted = 5;

                #Iterate through each GT column i.e. sample
                for my $i (9 .. $#{$record_ref}) {
                    my @gt_fields = split /:/, $record_ref->[$i];

                    #Update GT and add PGT:PID:PS
                    $gt_fields[$gt_i] =~ s/\//|/;
                    $gt_fields[$inserted] .= ":$gt_fields[$gt_i]:$pid:$ps";

                    $record_ref->[$i] = join ":", @gt_fields;
                    $record_ref->[$i] =~ s/:+/:/g;
                }
            }
        }
    ' $DONT_KEEP_NESTED_SBS |
    #Sort the records by coordinates
    perl -pe 's/^chrX/chr23/; s/^chrY/chr24/; s/^chrM/chr25/' |
    sort -k1,1V -k2,2V -k4,4V -k5,5V -T $TEMP_DIR |
    perl -pe 's/^chr23/chrX/; s/^chr24/chrY/; s/^chr25/chrM/' |
    #Add header
    cat <(bcftools view -h $INPUT_FILE | grep -vP "^##INFO=<ID=MBS_ID,") - |
    #Output as VCF.GZ
    bcftools view -Oz -o $OUTPUT_FILE

#Index the output VCF.GZ
tabix -p vcf $OUTPUT_FILE

