#!/usr/bin/perl -w
#
# A simple script for extracting pfam info.
#
# Author: Yuqian Jiang
# Created: 2023-06-14
# Version: 1.1.0
#
# Change logs:
# Update: 23-06-14: The initial version.
# Update: 23-06-28: Bug fixes: Col3 DESC will remove the first _, add a judge of the input file.
# Update: 23-06-28: Add GetOptions part and use MyFileIO module.
# Update: 24-01-05: unshift to print out the headline to the first line.

use strict;
use warnings;
use autodie;
use Path::Tiny;
use Getopt::Long;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

hmm_info.pl - extract info from .hmm files

=head1 SYNOPSIS

    perl hmm_info.pl -i <input_hmm> [options]
    A script for extracting info from pfam-hmm database.

      Options:
        --input    -i STR   input hmm file
        --out      -o STR   output file, default: STDOUT
        --help     -h       show help message

=cut

GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'in|i=s'    => \( my $input ),
    'out|o=s'   => \( my $out = 'stdout' ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($input) -> is_file ) {
    die "Error: can't find file [$input]";
}
elsif ( !$input =~ /\.hmm.*$/ ) {
    die ("Error: input a .hmm file please.");
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my ($readline, $line_num);
my @output = ();

my $fh_in = raid::MyFileIO::getInputFilehandle($input);
while( <$fh_in> ) {
    $line_num++;
    if ( $_ =~ /\/\// ) {
        if ( $readline =~ /NAME\s+?(\w.*?)\nACC\s+?(\w.+?)\nDESC\s+?(.+?)\nLENG\s+(\d+?)\n/ ){
            my $NAME = $1;
            my $ACC = $2;
            my $DESC = $3;
            my $LENG = $4;
            $DESC =~ s/\s/_/g;
            $DESC =~ s/^_//;
            my $for_print = "$NAME\t$ACC\t$DESC\t$LENG";
            push @output, $for_print;
            $readline = "";
        }
        else {
            print STDERR "Warning: HMM info extraction failed at line $line_num.\n";
            $readline = "";
        }
    }
    elsif ( $_ =~ /^\w/) {
        $readline .= $_;
    }
    else {
        next;
    }
}

close $fh_in;

my $head = "NAME\tACC\tDESC\tLENG";
unshift @output, $head;

raid::MyFileIO::print_out(\@output, $out);

__END__
