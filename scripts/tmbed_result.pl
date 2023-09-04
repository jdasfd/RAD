#!/usr/bin/perl -w
#
# A script aiming to converting TMbed result .pred to .bed format for better viewing.
#
# Author: Yuqian Jiang
# Created: 2023-07-07
# Version: 1.0.0
#
# Change logs:
# Version 1.0.0 23-07-08: The initial version.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use AlignDB::IntSpan;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

tmbed_result.pl - convert TMbed result .pred to .bed format

=head1 SYNOPSIS

    tmbed_result.pl - convert TMbed result .pred to .bed format

    perl tmbed_result.pl -f <fasta_file> [options]

    Options:
        --help          -h          brief help message
        --fasta         -f  STR     pep fasta file for quality control
        --out           -o  STR     output file name, default: STDOUT

=cut

GetOptions(
    'help|h'      => sub {Getopt::Long::HelpMessage(0)},
    'fasta|f=s'   => \(my $fasta_file),
    'out|o=s'     => \(my $out = 'stdout'),
) or Getopt::Long::HelpMessage(1);

if ( !defined $fasta_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($fasta_file) -> is_file ) {
    die "Error: can't find file [$fasta_file]";
}

#----------------------------------------------------------#
# main
#----------------------------------------------------------#

my %INDEXs;
my @bed_format;

raid::MyFileIO::read_pred(\%INDEXs, $fasta_file);

for my $id (keys %INDEXs) {
    my $sp_set = AlignDB::IntSpan -> new;
    my $tmd_o2i = AlignDB::IntSpan -> new;
    my $tmd_i2o = AlignDB::IntSpan -> new;

    my $index = $INDEXs{$id};
    my $length = length($index);

    for ( my $i = 0; $i < $length; $i++ )  {
        my $pos = $i+1;
        my $char = substr $index, $i, 1;
        if ( $char =~ /S/ ) {
            $sp_set -> add ($pos);
        }
        elsif ( $char =~  /[h|b]/ ) {
            $tmd_o2i -> add ($pos);
        }
        elsif ( $char =~ /[H|B]/ ) {
            $tmd_i2o -> add ($pos);
        }
    }

    unless ( $sp_set -> is_empty ) {
        my $run = $sp_set -> runlist;
        my @range = split/,/, $run;
        for (@range) {
            $_ =~ /(\d+)-(\d+)/;
            my $for_print = "$id\t$1\t$2\tSig_Pep";
            push @bed_format, $for_print;
        }
    }

    unless ( $tmd_o2i -> is_empty ) {
        my $run = $tmd_o2i -> runlist;
        my @range = split/,/, $run;
        for (@range) {
            $_ =~ /(\d+)-(\d+)/;
            my $for_print = "$id\t$1\t$2\tTMD_o2i";
            push @bed_format, $for_print;
        }
    }
    unless ( $tmd_i2o -> is_empty) {
        my $run = $tmd_i2o -> runlist;
        my @range = split/,/, $run;
        for (@range) {
            $_ =~ /(\d+)-(\d+)/;
            my $for_print = "$id\t$1\t$2\tTMD_i2o";
            push @bed_format, $for_print;
        }
    }
}

raid::MyFileIO::print_out(\@bed_format, $out);

__END__
