#!/usr/bin/perl -w
#
# MyFileIO.pm -- Process of files
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.0.0
#
# Change logs:
# Version 1.0.0 2023-06-08: Initial version

=head1 NAME

raid::MyFileIO - Local perl module for operating files

=head1 SYNOPSIS

    use raid::MyFileIO qw();

=cut

package raid::MyFileIO;
use strict;
use warnings;
use autodie;
use Bio::Seq;
use Bio::SeqIO;

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

=head2 read_pred

    About : Reading sequences from tmbed result (.pred)
    Usage : my @seq_id_list = read_pred(\%SEQs, $filename);
     Args : Reference to a hash to hold all the sequences;
            Filename (.pred format);
  Returns : Array of sequence IDs

=cut
sub read_pred {
    my ($rf_seq, $in) = @_;
    my @id_list;

    my $fh = getInputFilehandle($in);
    my $fas = do { local $/; <$fh> };

    my @fas = split /\>/, $fas;
    $fas[0] =~ s/^\>//;

    for my $str (@fas)
    {
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
            if ( ref($rf_seq) eq 'HASH' ) {
                $rf_seq -> {$id} = $tmindex;
            }
        }
        # save id
        push @id_list, $id;
    }
    return (@id_list);
}

1;

=head1 VERSION

1.0.0

=head1 AUTHOR

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

Special thanks to Nowind, noitulove9891@gmail.com

=cut
