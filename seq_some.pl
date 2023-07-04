#!/usr/bin/perl -w
#
# Inspired by faops
# A simple script for extracting sequences.
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

seq_some.pl - extract sequence by ID

=head1 SYNOPSIS

    perl seq_some.pl -f <fasta_file> -l <list_file>
      Options:
        --help          -h          brief help message
        --fasta         -f  STR     fasta file
        --list          -l  STR     ID list file

    perl seq_some.pl -f seq.fa -l seq_ID.lst

=cut

GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'fasta|f=s' => \( my $fasta_file ),
    'list|l=s'  => \( my $list_file ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file)->is_file ) {
    die "Error: can't find file [$fasta_file]";
}

if ( !defined $list_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($list_file)->is_file ) {
    die "Error: can't find file [$list_file]";
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

open my $l_in, '<', $list_file;

while (<$l_in>) {
    chomp;
    if ( exists $SEQUENCE{$_} ) {
        print ">$_\n";
        print "$SEQUENCE{$_}\n";
    }
    else {
        next;
    }
}

close $l_in;

__END__
