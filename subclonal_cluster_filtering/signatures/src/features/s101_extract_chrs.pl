#!/usr/bin/env perl

#Create sequence only files for each chromosome in given list for easier access
#later

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#Options for files to open (list of desired chrs, and genome build fasta), and output
#folder
    GetOptions(
        'list|l=s' => \my $list,
        'build|b=s' => \my $build,
        'output|o=s' => \my $out,
    ) or die "Invalid options";
    if (!$list || !$build || !$out) {
        print "Missing arguments";
        exit 0;
    }

#Read in list of desired chromosomes
    open(LIST, "<$list");
    my %chrs_desired;
    while (<LIST>) {
        chomp;
        $chrs_desired{$_} = 1;
    }
    close(LIST);

#Read in genome build file
    open(GENOME, "<$build");
    my $line;
    my $correct_chr = 0;

    while (<GENOME>) {
        chomp;

        #Check if start of new chromosome
        if ((substr $_, 0, 1) eq ">") {
            close;
            $correct_chr = 0;

            #Check if chromosome is desired
            (my $line) = $_ =~ />([^>\s]+)(\s|$)/;
            if (exists $chrs_desired{$line}) {
                open(OUTPUT, ">" . $out . $line);
                $correct_chr = 1;
            }
        } elsif ($correct_chr) {
            #Add to the file
            print OUTPUT $_;
        }
    }

    close;

