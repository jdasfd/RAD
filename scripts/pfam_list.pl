#!/usr/bin/perl -w
# Author: Yuqian Jiang
# A simple script for extracting pfam info.
# Date: latest update on 6/14/2023

use strict;
use warnings;
use autodie;
use Path::Tiny;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

pfam_list.pl - extract info from pfam.hmm database

=head1 SYNOPSIS

    perl pfam_list.pl <input_file>
    A script for extracting info from pfam-hmm database.

    Options:
    help: show help message
=cut

my $input = $ARGV[0];

die usage() if $input eq "help";
if ( !defined $input ) {
    die ("Input a file please.")
}
elsif {
    !path($input)->is_file;
    die ("Error: can't find file [$input]");
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my $readline;
my @output;

open my $fh_in, "<", $input;
while( <$fh_in> ) {
    if ( $_ =~ /\\\\/ ) {
        $readline =~ /NAME\s+?(\S.+?)\nACC\s+?(\S.+?)\nDESC\s+(\S.+?)\nLENG\s+(\d+?)\n/;
        my $NAME = $1;
        my $ACC = $2;
        my $DESC = $3;
        my $LENG = $4;
        $DESC = s/\s/_/g;
        my $for_print = "$NAME\t$ACC\t$DESC\t$LENG\n";
        push @output, $for_print;
        $readline = "";
    }
    else {
        $readline .= $_;
    }
}

END{
    print "NAME\tACC\tDESC\tLENG\n";
    foreach (@output) {
        print $_;
    }
}

sub usage{
    my $help_str = <<"EOF";
    perl pfam_list.pl <input_file>
    A script for extracting info from pfam-hmm database.

    Options:
    help: show help message
EOF
    return $help_str;
}
