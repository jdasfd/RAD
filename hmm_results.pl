#!/usr/bin/perl -w
#
# A simple script extract hmmscan txt output results.
#
# Author: Yuqian Jiang
# Created: 2023-06-28

# Change logs:
# Version: 1.0.0 23-06-14: The initial version, move from the RLK_family repo.
# Version: 1.0.1 23-06-29: Add more code for skipping empty query.

use strict;
use warnings;
use Path::Tiny;
use Bio::SearchIO;
use Getopt::Long;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

hmm_result.pl - transform hmm results to tsv format

=head1 SYNOPSIS

    perl hmm_result.pl -i <input_file> [stdout|-o output_file]
    A script for hmm results transforming.

Options:
    -i,--in STR     input file name (in txt format only)
    -o,--out STR    output file, default: STDOUT
    -h,--help       help information

=cut

GetOptions(
    "i|in=s"        => \(my $input),
    "o|out=s"       => \(my $out = 'stdout'),
    "h|help"        => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: please input a file.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($input) -> is_file ) {
    print STDERR "Error: can't find file [$input].\n";
    die Getopt::Long::HelpMessage(1);
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my (@for_print, @empty_domain);

my $searchio = Bio::SearchIO -> new (
                                        -format     => 'hmmer',
                                        -version    => 3,
                                        -file       => $input,
                                    );

while ( my $result = $searchio -> next_result() ) {
    while ( my $hit = $result -> next_hit ) {
            while ( my $hsp = $hit -> next_hsp ){
                    my $query = $result -> query_name;
                    my $hit_name = $hit -> name();
                    my $evalue = $hsp -> evalue();
                    my $start = $hsp -> start('query');
                    my $end = $hsp -> end('query');
                    my $for_print = "$query\t$hit_name\t$evalue\t$start\t$end";
                    push @for_print, $for_print;
            }
    }
}

raid::MyFileIO::print_out(\@for_print, $out);

__END__
