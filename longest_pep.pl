#!/usr/bin/perl -w
#
# longest_pep.pl -- extract longest protein sequences among pep
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.0.0
#
# Change logs:
# Version 1.0.0 06/08/23: The initial version.
# Version 1.0.1 23-06-21: Updated get_represent_trans sub.
# Version 1.0.2 23-06-21: Move get_represent_trans to perl module MyFileIO.
# Version 1.1.0 23-06-25: Realize longest transcripts function and move them into perl module MyFileIO.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use List::Utils qw(all);
use Bio::Seq;
use Bio::SeqIO;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

longest_pep.pl - extract longest peptide or protein files in fasta format via gff

=head1 SYNOPSIS

    perl longest_pep.pl -f <fasta_file> -g <gff_file> [options]
      Options:
        --fasta         -f  STR     input protein fasta file
        --gff           -g  STR     gff annotation file
        --help          -h          brief help message

    perl longest_pep.pl -f species_pep.fa -g species.gff

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'fasta|f=s'  => \( my $fasta_file ),
    'gff|g=s'    => \( my $gff_file ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    print STDERR "Error: please supply a protein fasta.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file) -> is_file ) {
    die "Error: can't find file [$fasta_file].";
}

if ( !defined $gff_file ) {
    print STDERR "Error: please supply an annotation gff3 format file.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($gff_file) -> is_file ) {
    die "Error: can't find file [$gff_file].";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my %LONGEST;
my $gff_in = raid::MyFileIO::getInputFilehandle($gff_file);

raid::MyFileIO::get_longest_trans(\%LONGEST, $gff_in);

print "Gene\tLongest_trans\n";

for my $gene_id (keys %LONGEST) {
    print "$gene_id\t$LONGEST{$gene_id}\n";
}
