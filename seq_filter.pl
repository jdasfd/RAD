#!/usr/bin/perl -w
#
# A script for discarding unqualified sequences.
#
# Author: Yuqian Jiang
# Created: 2023-07-05
# Version: 1.0.2
#
# Change logs:
# Version 1.0.0 23-07-05: The initial version.
# Version 1.0.1 23-07-05: Bug fixes. Remove arg: --changestop. Default change id and stop codon.
# Version 1.0.2 23-07-05: Add new module Array::Utils for intersect fuction and remove intersect sub.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use Bio::Seq;
use Bio::SeqIO;
use Array::Utils qw(:all);
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use rad::MyFileIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

seq_filter.pl - discard unqualified sequences of proteins (peptides)

=head1 SYNOPSIS

    seq_filter.pl - discard unqualified sequences of proteins (peptides)

    perl seq_filter.pl -f <fasta_file> [options]

    Notice: Stop codon will be replaced to * (default)
            ID with : will be replaced to _ (default)

    Options:
        --help          -h          brief help message
        --fasta         -f  STR     pep fasta file for quality control
        --longest       -l  NUM     remove seqs longer than this number (aa)
        --prestop                   remove seqs with premature stop codon (pep: * as stop codon)
        --out           -o  STR     output file name, default: STDOUT

=cut

GetOptions(
    'help|h'      => sub {Getopt::Long::HelpMessage(0)},
    'fasta|f=s'   => \(my $fasta_file),
    'longest|l=i' => \(my $length),
    'prestop'     => \(my $prestop),
    'out|o=s'     => \(my $out = 'stdout'),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file) -> is_file ) {
    die "Error: can't find file [$fasta_file]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

our %SEQUENCE;
my @filter_seq;
my @for_print;

# read fasta file
my $seqIOobj = Bio::SeqIO -> new(
    -file       =>  "$fasta_file",
    '-format'   =>  'Fasta'
);

while((my $seqobj = $seqIOobj -> next_seq())) {
    my $id = $seqobj -> id();
    my $seq = $seqobj -> seq();
    # modified ids and seqs here
    $id =~ s/:/_/g;
    $seq =~ s/\./\*/g;
    $SEQUENCE{$id} = $seq;
}

#----------------------------------------------------------#
# main
#----------------------------------------------------------#

# remove prestop & length
if ( defined $prestop && defined $length ) {
    my @fil_len = &filter_len(\%SEQUENCE, $length);
    my @fil_pre = &filter_pre_stop(\%SEQUENCE);
    @filter_seq = intersect(@fil_len, @fil_pre);
}
elsif ( defined $prestop ) {
    @filter_seq = &filter_pre_stop(\%SEQUENCE);
}
elsif ( defined $length ) {
    @filter_seq = &filter_len(\%SEQUENCE, $length);
}
else {
    @filter_seq = keys(%SEQUENCE);
}

for (@filter_seq) {
    my $print = ">$_\n$SEQUENCE{$_}";
    push @for_print, $print;
}

rad::MyFileIO::print_out(\@for_print, $out);

#----------------------------------------------------------#
# sub
#----------------------------------------------------------#

# Usage: @filter = filter_len(\%SEQUENCE, $length);
# return an array including all filtered seq ids
sub filter_len {
    my ($seq_hash, $len) = @_;
    my @filter;

    for my $id ( keys %{$seq_hash} ) {
        my $seq = $seq_hash -> {$id};
        my $seq_len = length($seq);
        if ( $seq_len <= $len ) {
            push @filter, $id;
        }
        else {
            next;
        }
    }

    return @filter;
}

# Usage: @filter = filter_pre_stop(\%SEQUENCE);
# return an array including seq ids discarding premature stop codon
sub filter_pre_stop {
    my ($seq_hash) = @_;
    my @filter;

    for my $id ( keys %{$seq_hash} ) {
        my $seq = $seq_hash -> {$id};
        $seq =~ s/\*$//;
        unless ( $seq =~ /\*/ ) {
            push @filter, $id;
        }
    }

    return @filter;
}

__END__
