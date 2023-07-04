#!/usr/bin/perl -w
#
# Inspired by faops
# A simple script for renaming all sequences.
#
# Author: Yuqian Jiang
# Created: 2023-03-27

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

seq_replace.pl - extract sequence by ID

=head1 SYNOPSIS

    perl seq_replace.pl -f <fasta_file> -t <tsv_file>
      Options:
        --help          -h          brief help message
        --fasta         -f  STR     fasta file
        --tsv           -t  STR     Replace name in tsv format (col1: raw name, col2: replacing name)

    perl seq_replace.pl -f seq.fa -t seq_ID.tsv

=cut

GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'fasta|f=s' => \( my $fasta_file ),
    'tsv|t=s'   => \( my $tsv_file ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file)->is_file ) {
    die "Error: can't find file [$fasta_file]";
}

if ( !defined $tsv_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($tsv_file)->is_file ) {
    die "Error: can't find file [$tsv_file]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

our %SEQUENCE;

# read fasta file
my $seqIOobj = Bio::SeqIO -> new(
                                -file       =>  "$fasta_file",
                                '-format'   =>  'Fasta'
                                );

while((my $seqobj = $seqIOobj -> next_seq())) {
    my $id = $seqobj -> id();
    my $seq = $seqobj -> seq();
    $SEQUENCE{$id} = $seq;
}

open my $t_in, '<', $tsv_file;

while (<$t_in>) {
    chomp;
    my @list = split/\t/, $_;
    if ( exists $SEQUENCE{$list[0]} ) {
        print ">$list[1]\n";
        print "$SEQUENCE{$list[0]}\n";
    }
    else {
        next;
    }
}

close $t_in;

__END__
