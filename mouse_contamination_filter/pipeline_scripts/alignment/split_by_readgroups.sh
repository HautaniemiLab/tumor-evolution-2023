#!/bin/bash

#Shell script for splitting BAM/SAM with multiple read groups to individual BAMs
#by read groups. Saves on I/O by reading input SAM/BAM only once.

#Required inputs are input SAM/BAM, output table and output BAM directory path.
#Optional input is custom sample key for output table.

#Requires SAMtools

USAGE="Usage: $(basename "$0") -i <input SAM/BAM> -o <output table> -d <output BAM dir> [OPTIONS]

Mandatory options:
    -i, --input         input BAM
    -d, --directory     output directory path for output BAMs
    -o, --output        output table of read groups and output BAM paths

Other options:
    -h, --help          show help
    -l, --link          use soft link if there is only one read group
    -s, --sample-key    sample name for output table
    -f, --include       only include reads with given SAM flags.
                        Disables linking
    -F, --exclude       remove reads with given SAM flags. Disables
                        linking
    -x, --do-not-index  do not index output BAM
    -p, --pipe-out      custom command for piping SAMs for each read
                        group. Useful albeit much slower for generating
                        sorted outputs. Defaults 'samtools view -b' to
                        output BAMs. Changing this disables linking"

#Parameters
    #Defaults
    SAMPLE=""
    LINK="FALSE"
    INCLUDE=""
    EXCLUDE=""
    INDEX="TRUE"
    PIPE_OUT="samtools view -b"

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h | --help )           >&2 echo "$USAGE"
                                    exit ;;
            -i | --input )          shift
                                    INPUT=$1 ;;
            -o | --output )         shift
                                    OUT_TBL=$1 ;;
            -d | --directory )      shift
                                    OUT_DIR=$1 ;;
            -l | --link )           LINK="TRUE" ;;
            -s | --sample-key )     shift
                                    SAMPLE=$1 ;;
            -f | --include )        shift
                                    INCLUDE="-f $1" ;;
            -F | --exclude )        shift
                                    EXCLUDE="-F $1" ;;
            -X | --do-not-index )   INDEX="FALSE" ;;
            -p | --pipe-out )       shift
                                    PIPE_OUT="$1" ;;
            * )                     >&2 echo "Illegal option: $1"
                                    >&2 echo "$USAGE"
                                    exit 1
        esac
        shift
    done

#Disable linking for specific options
if [[ $INCLUDE != "" || $EXCLUDE != "" || $PIPE_OUT != "samtools view -b" ]]; then
    LINK="FALSE"
fi

#Make sure obligatory arguments are given
if [[
    -z $INPUT ||
    -z $OUT_DIR ||
    -z $OUT_TBL
]]; then
    >&2 echo "Missing obligatory arguments!"
    >&2 echo "$USAGE"
    exit 1;
fi

#Create output directory if it doesn't exist
if [[ ! -d $OUT_DIR ]]; then
    mkdir $OUT_DIR
fi

set -euo pipefail

#Get sample name from first read group if not given
if [[ $SAMPLE == "" ]]; then
    SAMPLE=$(
        samtools view -H $INPUT |
            grep "^\@RG" | head -1 |
            perl -ne 'print $1 if m/SM:([^\t\n]+)/'
    )
fi

#Compute number of read groups
NUM_SAMPLES=$(samtools view -H $INPUT | grep "^\\@RG" | wc -l)

#Split BAM/SAM by read groups to respective files, reading ONCE
if [[ $NUM_SAMPLES -gt 1 || $LINK != "TRUE" ]]; then
    samtools view -h $INPUT $INCLUDE $EXCLUDE |
        perl -e '
            use warnings; use strict;

            my ($sample_key, $out_table, $out_folder, $pipe_out) = @ARGV;

            my @pre_rg_lines;
            my @rg_lines;
            my @rg_fields;
            my %rg_fields;
            my %rg_files;

            my $phase = 0;

            while (<STDIN>) {
                #Store header lines preceding @RG
                if ($phase == 0) {
                    if (not m/^\@RG/) {
                        push @pre_rg_lines, $_;
                    } else {
                        $phase++;
                    }
                }

                #Create read group file for each @RG line and create respective line to table
                if ($phase == 1) {
                    if (m/^\@RG/) {
                        #Add RG line
                        push @rg_lines, $_;

                        chomp;

                        #Add RG fields
                        for my $rg_col (split /\t/) {
                            if ($rg_col =~ m/^([A-Z]{2}):/) {
                                if (not exists $rg_fields{$1}) {
                                    push @rg_fields, $1;
                                }

                                $rg_fields{$1} = 1;
                            }
                        }
                    } else {
                        #Create output table and print header
                        open TABLE, ">", $out_table or die "$!";
                        print TABLE "sample\t" . join("\t", @rg_fields) . "\tfile\n";

                        #Add metadata line to table and create file handle for each read group
                        for my $rg_line (@rg_lines) {
                            #Get RG ID
                            $rg_line =~ /ID:([^\t\n]+)/;
                            my $id = $1;

                            my $file_path = "$out_folder/$id.bam";

                            #Add read group file to table
                            print TABLE "$sample_key";

                            #Print RG field values
                            for my $rg_field (@rg_fields) {
                                if ($rg_line =~ m/${rg_field}:([^\t\n]+)/) {
                                    print TABLE "\t$1";
                                } else {
                                    print TABLE "\tNA";
                                }
                            }
                            print TABLE "\t$file_path\n";

                            #Create output file handle with piping to samtools BAM conversion
                            open $rg_files{$id}, "| $pipe_out > $file_path" or die "$!";

                            #Write header to new BAM
                            for my $line (@pre_rg_lines) {
                                print {$rg_files{$id}} $line;
                            }
                            print {$rg_files{$id}} $rg_line;
                        }

                        $phase++;
                    }
                }

                #Write rest of the header to each read group file
                if ($phase == 2) {
                    if (m/^\@/) {
                        for my $id (keys %rg_files) {
                            print {$rg_files{$id}} $_;
                        }
                    } else {
                        $phase++;
                    }
                }

                #Write reads to respective read group file
                if ($phase == 3) {
                    / RG:Z:([^\t]+) /x;
                    print {$rg_files{$1}} $_;
                }
            }

            close TABLE;
            for my $id (keys %rg_files) {
                close $rg_files{$id};
            }
        ' $SAMPLE $OUT_TBL $OUT_DIR "$PIPE_OUT"

    #Index output BAMs
    if [[ $INDEX == "TRUE" ]]; then
        while read FILE; do
            samtools index $FILE &
        done < <(perl -ane 'if ($. > 1) {print "$F[$#F]\n"}' $OUT_TBL)

        wait
    fi
else
    #Create soft links
    FILE_PATH=$OUT_DIR/$(basename $INPUT)
    ln -s $INPUT $FILE_PATH
    if [[ $INDEX == "TRUE" && -f $INPUT.bai || -f $(dirname $INPUT)/$(basename $INPUT .bam).bai ]]; then
        ln -s $(dirname $INPUT)/$(basename $INPUT .bam)*.bai $OUT_DIR/
    fi

    #Create output table
    samtools view -H $INPUT |
        perl -e '
            use warnings; use strict;

            my ($sample_key, $out_table, $file_path) = @ARGV;

            #Look for RG line
            while (<STDIN>) {
                if (m/^\@RG/) {
                    #Create output table and print header
                    open TABLE, ">", $out_table or die "$!";
                    print TABLE "sample";

                    my @rg_vals;

                    for my $rg_col (split /\t/) {
                        if ($rg_col =~ m/^([A-Z]{2}):([^\t\n]+)/) {
                            print TABLE "\t$1";

                            push @rg_vals, $2;
                        }
                    }
                    print TABLE "\tfile\n";

                    print TABLE "$sample_key\t" . join("\t", @rg_vals) . "\t$file_path\n";

                    close TABLE;

                    last;
                }
            }
        ' $SAMPLE $OUT_TBL $FILE_PATH
fi

