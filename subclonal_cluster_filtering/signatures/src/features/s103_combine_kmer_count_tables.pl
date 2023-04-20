#!/usr/bin/env perl

#Combines kmer counts tables from individual tables (i.e. chromosomes). The
#tables to be combined are given in a list file. Combines by summing respective
#fields using hash.
#Outputs to stdout

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#Requires the list of files
    GetOptions(
        'list|l=s' => \my $list,
    ) or die "Invalid options";
    if (!$list) {
        print "Missing arguments";
        exit 0;
    }

#Variables; hash for kmer-count, respective variables per line as well as first line
#flag and first line
    my %kmer_counts;
    my $kmer;
    my $count;
    my $headered = 0;
    my $first_line;

#Open the list file and go through the tables in the list, adding up the kmers
#encountered
    open LIST, "<$list";

    while (<LIST>) {
        chomp;

        #Open the file
        open TBL, ("<" . $_);

        #Handle the header
        $first_line = <TBL>;
        if (!$headered) {
            print $first_line;
            $headered = 1;
        }

        #Iterate through data rows
        while (<TBL>) {
            chomp;

            #Get kmer and count from line
            ($kmer, $count) = split "\t";
            #Add value to hash
            if (exists($kmer_counts{$kmer})) {
                $kmer_counts{$kmer} += $count;
            } else {
                $kmer_counts{$kmer} = $count;
            }
        }

        close TBL;
    }

    close LIST;

#Print the results, i.e. kmer-count pairs
    for my $key (sort keys %kmer_counts) {
        print "$key\t$kmer_counts{$key}\n";
    }

