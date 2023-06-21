#!/usr/bin/perl -w
#
#   longest_pep.pl -- extract longest protein sequences among pep
#
#   Author: Yuqian Jiang
#   Created: 2023-06-08
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 06/08/23: The initial version.
#   Version 1.0.1 06/21/23: Updated get_represent_trans sub.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use Bio::Seq;
use Bio::SeqIO;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;
use Data::Dumper;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

longest_pep.pl - extract longest peptide or protein files in fasta format via gff

=head1 SYNOPSIS

    perl longest_pep.pl -f <fasta_file> -g <gff_file>
      Options:
        --help          -h          brief help message
        --input         -i  STR     input protein fasta file
        --gff           -g  STR     gff annotation file

    perl longest_pep.pl -i species_pep.fa -g species.gff

=cut

GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
#    'input|i=s' => \( my $input_file ),
    'gff|g=s'   => \( my $gff_file ),
) or Getopt::Long::HelpMessage(1);

=pod
if ( !defined $input_file ) {
    print STDERR "Error: please supply a protein fasta.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($input_file)->is_file ) {
    die "Error: can't find file [$input_file].";
}
=cut

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
