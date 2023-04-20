#!/usr/bin/env perl

#Compute the total number of kmers within the genome, using given list of
#chromosomes to count the kmers from. Possibility to force pyrimidine on certain
#position or flip dinucleotides as in Alexandrov et al. 2020.
#Outputs to stdout

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#Requires either file name of chrs list and directory with chrs (including
#forward slash '/') or chromosome file.
#Options include position with forced pyrimidine (default off with non-positive
#values), dinucleotide flip following Alexandrov et al. 2020, and kmer size.
#Intervals BED file can be given to count only from subsets of whole
#chromosomes, and padding can be supplied to extend the ends of intervals.
    my $pyr = -1;
    my $kmer_size = 3;
    my $padding = 0;
    GetOptions(
        'list|l=s' => \my $list,
        'dirchr|d=s' => \my $chr_dir,
        'chr|c=s' => \my $chr_file,
        'pyrimidine|p|=i' => \$pyr,
        'dinuc_flip|f' => \my $dinuc_flip,
        'kmer_size|k=i' => \$kmer_size,
        'intervals=s' => \my $intervals_file,
        'padding=i' => \$padding,
    ) or die "Invalid options";
    if (!($list and $chr_dir) and !$chr_file) {
        print "Missing arguments!\n";
        exit 0;
    }
    if ($kmer_size < 1) {
        print "kmer size must be at least 1!\n";
        exit 0;
    }
    if ($pyr > $kmer_size) {
        print "Forced pyrimide position must be within [1, kmer_size], smaller values imply off\n!";
        exit 0;
    }
    if ($kmer_size != 2 and $dinuc_flip) {
        print "Dinucleotide flipping following Alexandrov et al. 2018 only available for 2-mers\n!";
        exit 0;
    }

#Variables; hash containing kmer count and strings
    my $chr_seq_ref;
    my @interval;
    my %kmer_counts;
    my $kmer;

#Read specified chromosome file, if available, otherwise use the given list
    if ($chr_file) {
        #Count kmers in chromosome
        $chr_seq_ref = read_chromosome($chr_file);
        count_kmers($chr_seq_ref, \%kmer_counts);
    } else {
        #Process chromosome list
        my %chromosomes = process_chromosome_list($list, $chr_dir);

        #Count by intervals if such file has been given, otherwise simply interate through all chromosomes
        if ($intervals_file) {
            my $curr_chr = "";

            open INTERVALS, "<$intervals_file";

            while (<INTERVALS>) {
                chomp;
                @interval = split /\t/;

                #Read in new chromosome if encountered
                if ($curr_chr ne $interval[0]) {
                    $curr_chr = $interval[0];
                    $chr_seq_ref = read_chromosome($chromosomes{$curr_chr});
                }

                #Count kmers in the interval, convert BED end position to inclusive position
                count_kmers($chr_seq_ref, \%kmer_counts, $interval[1], $interval[2]-1);
            }

            close INTERVALS;
        } else {
            #Count kmers in each chromosome
            foreach my $chr_file (values %chromosomes) {
                $chr_seq_ref = read_chromosome($chr_file);
                count_kmers($chr_seq_ref, \%kmer_counts);
            }
        }
    }

#Convert to pyrimidine representation if flagged
    if ($pyr > 0) {
        for my $key (keys %kmer_counts) {
            #Check if the middle base is C or T
            if ((substr $key, $pyr-1, 1) !~ /[CT]/) {
                #supposed representation
                ($kmer = reverse $key) =~ tr/ACGT/TGCA/;
                #add the counts together
                $kmer_counts{$kmer} += $kmer_counts{$key};
                #remove undesired representation
                delete $kmer_counts{$key};
            }
        }
    }

#Convert doublets according according to the classification in Alexandrov et al. 2018:
# - REF 5' preference for non-symmetric doublets
#   - T > C > A > G
    if ($dinuc_flip) {
        #REF 5' preference hash table
        my %ref_5_pref = (T => 1, C => 2, A => 3, G => 4);
        my @bases;

        #Iterate through dinucleotides
        for my $key (keys %kmer_counts) {
            @bases = split "", $key;
            #Check if reverse complement is to be reported
            if ($ref_5_pref{$bases[1] =~ tr/ACGT/TGCA/r} < $ref_5_pref{$bases[0]}) {
                #supposed representation
                ($kmer = reverse $key) =~ tr/ACGT/TGCA/;
                #add the counts together
                $kmer_counts{$kmer} += $kmer_counts{$key};
                #remove undesired representation
                delete $kmer_counts{$key};
            }
        }
    }

#Print the results
    #header
    print "kmer\tcount\n";

    #kmer-count pairs
    for my $key (sort keys %kmer_counts) {
        print "$key\t$kmer_counts{$key}\n";
    }

#Subroutine for opening chromosome file and reading in the sequence
sub read_chromosome {
    #Arguments: chromosome file path
    my $chr_file = shift;

    #Read chrosome file
    open CHR, ("<" . $chr_file);
    my $chr_seq = <CHR>;
    close CHR;

    return \$chr_seq;
}

#Subroutine for going through chromosome list
sub process_chromosome_list {
    #Arguments: chromosome list file path
    (my $list, $chr_dir) = @_;

    my %chromosomes;

    #Read through the list
    open LIST, "<$list";

    while (<LIST>) {
        chomp;
        $chromosomes{$_} = $chr_dir . $_;
    }

    close LIST;

    return %chromosomes;
}

#Subroutine for counting kmers
sub count_kmers {
    #Arguments: references to (sequence string, kmer_count hash); 0-based start and end (inclusive) positions
    (my $chr_seq_ref, my $kmer_counts_ref, my $start, my $end) = @_;

    #If given interval, specify proper start and end positions with potential padding
    $start -= $padding if ($start);
    $end += 1 - $kmer_size + $padding if ($end);
    #Make sure the interval is valid at contig ends
    if (! $start or $start < 0) {
        $start = 0;
    }
    if (! $end or $end > length($$chr_seq_ref) - $kmer_size) {
        $end = length($$chr_seq_ref) - $kmer_size;
    }

    #Iterate through file and add frequencies to kmer counts
    for my $pos ($start .. $end) {
        $kmer = uc substr $$chr_seq_ref, $pos, $kmer_size;

        #Check if kmer is eligible
        if ($kmer !~ /[^ACGT]+/) {
            #Check if kmer is already in hash
            if (exists($kmer_counts_ref->{$kmer})) {
                $kmer_counts_ref->{$kmer}++;
            } else {
                $kmer_counts_ref->{$kmer} = 1;
            }
        }
    }
}

