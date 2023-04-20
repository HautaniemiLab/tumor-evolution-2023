#!/usr/bin/env perl

#Perl script to determine ID typing (strand) according to the classification in
#Alexandrov et al. ? (bioRxiv 2018):
# - INS / DEL
# - size 1, 2, 3, 4, 5+
# - repeat size 0 (+1), 1, 2, 3, 4, 5+
# - deletions with repeat size <2 checked for microhomology
#   - deletion size 2, 3, 4, 5+
#   - max. microhomology of either 5' or 3' directions, upto deletion size -1
#     homology size 1, (2), (3), (4), (5+)

use strict;
use warnings;

#Loop over input file (VCF without headers) given as stdin
    my @line;
    my $id_seq;
    my $start_i;
    my $type;
    my $id_size;
    my $count;

    while (<STDIN>) {
        #Read line
        chomp;
        @line = split "\t";

        #Get inserted/deleted sequence
        if ($line[5] eq "INS") {
            $id_seq = substr $line[3], 1, $line[6];
        } else {
            $id_seq = substr $line[2], 1, $line[6];
        }
        #Repeat start position
        $start_i = $line[6] - 1;

        #ID type to be concatenated
        $type = $line[5];

        #ID size for type concatenation
        $id_size = $line[6];
        if ($id_size >= 5) {$id_size = "5+"}

        #Check repeat count
        $count = repeat_count($id_seq, $line[7], $start_i, $line[6]);

        #Check length if mononucleotide indel
        if ($line[6] == 1) {
            #C(:G) or T(:A) base pair
            if ($id_seq =~ m/[CG]/) {
                $type .= "_C_1";
            } else {
                $type .= "_T_1";
            }
        } elsif ($line[5] eq "DEL" and $count == 1) {
            #Check if micro-homology
            #3' loop
            my $i;
            for ($i=0; $i<$line[6]-1; $i++) {
                #End at first mismatch
                if (substr($line[7], $start_i + $i, 1) ne substr($line[7], $start_i + $line[6] + $i, 1)) {
                    last;
                }
            }
            my $mh_size = $i;
            #5' loop
            for ($i=1; $i<$line[6]; $i++) {
                #End at first mismatch
                if (substr($line[7], $start_i - $i, 1) ne substr($line[7], $start_i + $line[6] - $i, 1)) {
                    last;
                }
            }
            if ($i - 1 > $mh_size) {$mh_size = $i - 1}
            #Save type and check if there was microhomology
            if ($mh_size > 0) {
                if ($mh_size >= 5) {$mh_size = "5+"}
                $type .= "_MH_${id_size}_$mh_size";
                $count = -1;
            } else {
                $type .= "_repeats_${id_size}";
            }
        } else {
            #Other INS/DEL
            $type .= "_repeats_${id_size}";
        }

        #Add repeat count for those without
        if ($count >= 0) {
            #Limit substitution; insertions
            if ($line[5] eq "INS") {
                if ($count >= 5) {$count = "5+"}
            } else {
                #Deletions, consider the deleted sequence as repeat number one?
                #$count--;
                #if ($count >= 5) {$count = "5+"}
                if ($count >= 6) {$count = "6+"}
            }
            #Concatenation
            $type .= "_$count";
        }

        #Print line
        print(join("\t", (@line[0 .. 4], $type, @line[5 .. $#line])) . "\n");
    }

#Subroutine for checking repeat count
sub repeat_count {
    #Arguments: id_seq, full sequence, start index, id length
    my ($id_seq, $full_seq, $start_i, $id_len) = @_;

    #For loop for counting repeats
    my $count;
    for ($count=0; $count<6; $count++) {
        #End counting once there is a mismatch
        if (substr($full_seq, $start_i + $count * $id_len, $id_len) ne $id_seq) {
            last;
        }
    }

    return $count;
}

