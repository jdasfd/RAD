#!/usr/bin/perl -w
#
# A simple script for translating cds into pep.
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
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;
use raid::OptSeq;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

seq_trans.pl - translate cds into pep

=head1 SYNOPSIS

    seq_trans.pl - translate cds into pep

    perl seq_trans.pl -i <cds_file>

    Options:
        --help          -h          brief help message
        --in            -i  STR     input cds file (fasta format)
        --out           -o  STR     output, default: STDOUT

=cut

GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'in|i=s'    => \( my $input_cds ),
    'out|o=s'   => \( my $out = 'stdout' ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $input_cds ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($input_cds)->is_file ) {
    die "Error: can't find file [$input_cds]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my @for_print;

# read fasta file
my $seqIOobj = Bio::SeqIO -> new(
    -file       =>  "$input_cds",
    '-format'   =>  'Fasta'
);

# print subseq
while ((my $seqobj = $seqIOobj -> next_seq())) {
    my $id = $seqobj -> id();
    my $seq = $seqobj -> seq();
    next unless ( length($seq)%3 == 0 );
    my $pep_seq = raid::OptSeq::codon_translate($seq);
    push @for_print, ">$id\n$pep_seq";
}

raid::MyFileIO::print_out(\@for_print, $out);

__END__
