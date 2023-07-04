#!/usr/bin/perl -w
#
# Inspired by faops
# A simple script for acquiring sub-domain sequence by truncating.
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

seq_trunc.pl - extract truncated sequence by giving a region

=head1 SYNOPSIS

    perl seq_trunc.pl -f <fasta_file> -t <tsv_file>
      Options:
        --help          -h          brief help message
        --fasta         -f  STR     fasta file
        --tsv           -t  STR     tsv file (col1: gene, col2: start, col3: end, col4: domain name)

    perl seq_trunc.pl -f seq.fa -l trunc_info.tsv

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

our %INFO;

open my $t_in, '<', $tsv_file;

while (<$t_in>) {
    chomp;
    my @array = split/\t/, $_;
    for ( my $i = 1; $i < 4; $i++ ) {
        push @{$INFO{$array[0]}}, $array[$i];
    }
}

close $t_in;

# read fasta file
my $seqIOobj = Bio::SeqIO -> new(
    -file       =>  "$fasta_file",
    '-format'   =>  'Fasta'
);

# print subseq
while ((my $seqobj = $seqIOobj -> next_seq())) {
    my $id = $seqobj -> id();
    if ( exists $INFO{$id} ) {
        my $dom_count = @{$INFO{$id}};
        my $dom_num = $dom_count/3;
        for ( my $i = 0; $i < $dom_num; $i++ ) {
            my $arr_num = $i*3;
            my $subseq = $seqobj -> subseq($INFO{$id} -> [$arr_num], $INFO{$id} -> [$arr_num+1]);
            my $domain = $INFO{$id} -> [$arr_num+2];
            my $domain_name = "$domain"."_$i";
            print ">$id","_$domain_name\n";
            print "$subseq\n";
        }
    }
    else {
        next;
    }
}

__END__
