#!/usr/bin/env perl

#Perl script to determine DBS typing (strand) according to the classification in
#Alexandrov et al. ? (bioRxiv 2018):
# - REF 5' preference for non-symmetric doublets
#   - T > C > A > G
# - ALT 5' preference for symmetric doublets, different for each REF 5'
#   - A: C > G > T
#   - C: T > G > A
#   - G: A > C > T
#   - T: G > C > A

use strict;
use warnings;

#Global variables
    #REF 5' preference hash table
    my %ref_5_pref = (T => 1, C => 2, A => 3, G => 4);
    #ALT 5' preference hash table
    my %alt_5_pref = (
        A => {C => 1, G => 2, T => 3},
        C => {T => 1, G => 2, A => 3},
        G => {A => 1, C => 2, T => 3},
        T => {G => 1, C => 2, A => 3}
    );

#Loop over input file (VCF without headers) given as stdin
    my @line;
    my @ref;
    my @alt;
    my $type;
    my $strand;

    while (<STDIN>) {
        #Read line
        chomp;
        @line = split "\t";

        #Decompose REF and ALT
        @ref = split "", $line[2];
        @alt = split "", $line[3];

        #Check if reverse strand needs to be reported
        if ($ref_5_pref{$ref[1] =~ tr/ACGT/TGCA/r} < $ref_5_pref{$ref[0]}) {
            #Reverse strand REF 5' is preferable
            $type = "$ref[1]$ref[0]>$alt[1]$alt[0]" =~ tr/ACGT/TGCA/r;
            $strand = "-";
        } elsif (($ref[1] =~ tr/ACGT/TGCA/r) eq $ref[0] && $alt_5_pref{$ref[0]}{$alt[1] =~ tr/ACGT/TGCA/r} < $alt_5_pref{$ref[0]}{$alt[0]}) {
            #Reverse strand ALT 5' is preferable when REF is an inverted repeat
            $type = "$ref[1]$ref[0]>$alt[1]$alt[0]" =~ tr/ACGT/TGCA/r;
            $strand = "-";
        } else {
            #Otherwise forward strand is preferred
            $type = "$ref[0]$ref[1]>$alt[0]$alt[1]";
            $strand = "+";
        }

        #Print line
        print(join("\t", (@line[0 .. 4], $type, $strand, @line[5 .. $#line])) . "\n");
    }

