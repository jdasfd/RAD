#!/usr/bin/env perl
#
# MyFileIO.pm -- Process of files
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.2.1
#
# Change logs:
# Version 1.0.0 2023-06-08: Initial version. Add function getInputFilehandle, read_pred.
# Version 1.1.0 2023-06-21: Add function get_represent_trans and rename it to get_longest_trans.
# Version 1.1.1 2023-06-25: Modified get_longest_trans codes.
# Version 1.1.2 2023-06-26: Fixes bugs: skip empty lines and accept transcript annotation.
# Version 1.2.0 2023-06-26: Add function print_out.
# Version 1.2.1 2023-07-08: Modified read_pred, return nothing now. And skip the first empty line.
# Version 1.3.0 2023-08-10: Add function read_hmm.

=head1 NAME

raid::MyFileIO - Local perl module for operating files

=head1 SYNOPSIS

    use raid::MyFileIO;

=cut

package raid::MyFileIO;

use strict;
use warnings;
use Bio::SearchIO;

=head1 METHODS

=head2 getInputFilehandle

      About : Open and return filehandles
      Usage : my $fh = getInputFilehandle($filename);
       Args : Filename
    Returns : Filehandle to the opened file

=cut
sub getInputFilehandle
{
    my ($in) = shift;

    my $expr = "-";
    ## read from STDIN
    if ( !defined $in || $in eq "-") {
        $expr = "-";
    }
    ## read from a tar gzip ball
    elsif ($in =~ /\.tar\.gz$/) {
        $expr = "tar -zxf $in -O |";
    }
    ## read from a tar bzip2 ball
    elsif ($in =~ /\.tar\.bz2$/) {
        $expr = "tar -jxf $in -O |";
    }
    ## read from a gzip file
    elsif ($in =~ /\.gz$/) {
        $expr = "gzip -dc $in |";
    }
    ## read from a bzip2 file
    elsif ($in =~ /\.bz2$/) {
        $expr = "bzip2 -dc $in |";
    }
    ## read from a zip file
    elsif ($in =~ /\.zip$/) {
        $expr = "unzip -p $in |";
    }
    else {
        $expr = "< $in";
    }

    open (my $fh, $expr) || die $!;

    return $fh;
}

=head2 print_out

    About : Output lines saved in the print array
    Usage : print_out(\@for_print, $out_file);
     Args : Reference to an array to hold all the content;
            Output file (stdout or filename);
  Returns : None

=cut
sub print_out{
    my ($print_array, $out) = @_;

    my $out_fh;

    if ( lc($out) eq "stdout" ) {
        $out_fh = *STDOUT;
    }
    else {
        open $out_fh, ">", $out;
    }

    for (@{$print_array}) {
        print {$out_fh} $_ . "\n";
    }

    close $out_fh;
}

=head2 read_pred

    About : Reading sequences from tmbed result (.pred)
    Usage : read_pred(\%INDEXs, $filename);
     Args : Reference to a hash to hold all the indexes with id;
            Filename (.pred format);
  Returns : Nothing

=cut
sub read_pred {
    my ($seq_hash, $in) = @_;

    my $fh = getInputFilehandle($in);
    my $file = do { local $/; <$fh> };

    my @fas = split /\>/, $file;
    for (@fas) {
        if ( $_ eq "" ) {
            shift @fas;
        }
    }

    for my $str (@fas) {
        $str =~ /^(.+?)\n/;
        my $id = $1;

        $str =~ s/^.+?\n//;
        my $line_num = $str =~ tr/\n/\n/;
        $line_num = $line_num/2;

        my ($seq, $tmindex);

        for (my $i = 0; $i < $line_num; $i++) {
            $str =~ /^(.+?)\n/;
            if ( !defined($seq) ) {
                $seq = $1;
            }
            else {
                $seq = $seq.$1;
            }
            $str =~ s/^.+?\n//;
        }

        $str =~ s/\n//g;
        $tmindex = $str;

        # length
        my $len_seq = length($seq);
        my $len_tm = length($tmindex);

        unless ( $len_seq != $len_tm ) {
            if ( ref($seq_hash) eq 'HASH' ) {
                $seq_hash -> {$id} = $tmindex;
            }
        }
    }
}

=head2 read_hmm

    About : Reading hmmscan result to a hash
    Usage : my @for_print = read_hmm($in);
     Args : Hmmscan result txt filename (.txt format);
  Returns : @for_print;

=cut
sub read_hmm {
    my ($in) = @_;
    my @for_print;

    my $searchio = Bio::SearchIO -> new (
        -format     => 'hmmer',
        -version    => 3,
        -file       => $in,
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

    return @for_print;
}

=head2 get_longest_trans

    About : Get the longest transcript as the representative transcript (NCBI format)
    Usage : get_longest_trans(\%LONGEST_TRANS, $gff_file);
     Args : Reference to a hash to contain all longest trans;
            gff filename (.gff3 format);
  Returns : longest transcripts hash [key: transcript id -> value: mRNA parent id]

=cut
sub get_longest_trans {
    my ($seq_hash, $in) = @_;
    my (%mRNAs_all, %CDS);
    my $fh = getInputFilehandle($in);
    while( <$fh> ) {
        next if (/^\#/ || /^\s+$/ || /^$/);
        chomp;
        my ($feature, $start, $end) = (split /\t/)[2,3,4];

        if ( $feature =~ /(mRNA|pseudogenic\_transcript|transcript)/ ) {
            m{
                ^(.*?)\s+.*?                   # Chromosome ID
                \s+(\d+)                       # Start Position
                \s+(\d+)\s+.*                  # End Position
                \s+(\-|\+)\s+\.                # Strand
                \s+ID\=(.*?)\;.*               # ID
                (Parent|Locus_id)\=(.*?)(;|$)  # Parent
            }x;

            my ($id, $parent) = ($5, $7);

            my $length = $end - $start + 1;

            if ( $mRNAs_all{$parent} -> {represent} ) {
                if ($mRNAs_all{$parent} -> {length} < $length) {
                    $mRNAs_all{$parent} -> {represent} = $id;
                    $mRNAs_all{$parent} -> {length} = $length;
                }
            }
            else {
                $mRNAs_all{$parent} -> {represent} = $id;
                $mRNAs_all{$parent} -> {length} = $length;
            }
        }
        elsif ( $feature =~ /CDS/ ) {
            m{
                \s+ID\=(.*?)\;.*
                (Parent|Locus_id)\=(.*?)(;|$)
            }x;

            my ($id, $parent) = ($5, $7);
            $CDS{$parent} = $id;
        }
    }

    for my $mRNA (keys %mRNAs_all)
    {
        my $represent_id = $mRNAs_all{$mRNA} -> {represent};

        $seq_hash -> {$represent_id} = $mRNA;
    }
}

1;

=head1 VERSION

1.3.0

=head1 AUTHOR

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

Special thanks to:
Qiang Wang, wangq@outlook.com
Nowind, noitulove9891@gmail.com

=cut
