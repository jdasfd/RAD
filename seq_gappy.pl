#!/usr/bin/perl -w
#
# A simple script for counting gaps in the sequence.
#
# Author: Yuqian Jiang
# Created: 2023-10-11
# Modified: 2023-12-20
# Decided to put it into RAID tool-box for filtering gappy sequences

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Path::Tiny;
use Bio::Seq;
use Bio::SeqIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

seq_gappy.pl - remove those gappy sequences

=head1 SYNOPSIS

    seq_gappy.pl - remove those gappy sequences

    perl seq_gappy.pl -f <fasta_file> -g <cut-off>

    Options:
        --help          -h          brief help message
        --fasta         -f  STR     fasta file
        --greater       -g  INT     gappy cut-off, filtering gap% greater than the value

=cut

GetOptions(
    'help|h'        => sub { Getopt::Long::HelpMessage(0) },
    'fasta|f=s'     => \( my $fasta_file ),
    'greater|g=s'   => \( my $cutoff ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file)->is_file ) {
    die "Error: can't find file [$fasta_file]";
}

if ( !defined $cutoff ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( $cutoff > 100 || $cutoff < 0) {
    die "Error: cutoff ranged from 0-100!";
}

# read fasta file
my $seqIOobj = Bio::SeqIO -> new(
    -file       =>  "$fasta_file",
    '-format'   =>  'Fasta'
);

# filter all sequence with gaps
while((my $seqobj = $seqIOobj -> next_seq())) {
    my $id = $seqobj -> id();
    my $seq = $seqobj -> seq();
    my $seq_len = length($seq);
    my $gap_len = $seq =~ tr/-/-/;
    my $ratio = $gap_len/$seq_len*100;
    $ratio = sprintf "%.3f", $ratio;
    if ( $ratio <= $cutoff ) {
        print ">$id\n";
        print "$seq\n";
    }
}

__END__
