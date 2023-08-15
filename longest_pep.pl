#!/usr/bin/perl -w
#
# longest_pep.pl -- extract longest protein sequences among pep
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.2.0
#
# Change logs:
# Version 1.0.0 23-06-08: The initial version.
# Version 1.0.1 23-06-21: Updated get_represent_trans sub.
# Version 1.0.2 23-06-21: Move get_represent_trans to perl module MyFileIO.
# Version 1.1.0 23-06-25: Realize longest transcripts function and move them into perl module MyFileIO.
# Version 1.2.0 23-06-26: Complete longest pep extraction part in two situation that without multi-transcripts.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use List::Util qw(all);
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
        --output        -o  STR     output file name, default: STDOUT
        --help          -h          brief help message

    perl longest_pep.pl -f species_pep.fa -g species.gff

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'fasta|f=s'  => \( my $fasta_file ),
    'gff|g=s'    => \( my $gff_file ),
    'output|o=s' => \( my $output = 'stdout' )
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

my (%LONGEST, %SEQUENCE);
my (@GFFLINES, @f_array, @for_print);

# This part code are judging different gff types. Many gff files do not contain multiple transcripts.
# We should judge the situation first, and then get the longest pep or transcript.

print STDERR "==> Start to extract longest pep.\n";

my $gff_in = raid::MyFileIO::getInputFilehandle($gff_file);

while ( <$gff_in> ) {
    chomp;
    next if (/^\#/ || /^\s+$/ || /^$/);
    push @GFFLINES, $_;
    my $features = (split /\t/)[2];
    push @f_array, $features;
}
close $gff_in;

my @f_uniq = List::Util::uniq @f_array;
if ( List::Util::any { $_ eq "mRNA" || $_ eq "pseudogenic_transcript" || $_ eq "transcript" } @f_uniq ) {
    if ( List::Util::any { $_ =~ /.*mRNA.*Parent.*/ } @GFFLINES ) {
        #raid::MyFileIO::get_longest_trans(\%LONGEST, $gff_file);
        print STDERR "==> This part is waiting to update QAQ!\n";
    }
    else {
        print STDERR "==> mRNA without [Parent;] note.\n";
        my @ID_list;
        for (@GFFLINES) {
            m{
                ^(.*?)\s+.*?                   # Chromosome ID
                \s+(\d+)                       # Start Position
                \s+(\d+)\s+.*                  # End Position
                \s+(\-|\+)\s+\.                # Strand
                \s+ID\=(.*?)\;.*               # ID
            }x;
            my $id = $5;
            push @ID_list, $id;
        }
        my $mRNA_num = @ID_list;
        print STDERR "==> Total $mRNA_num mRNA in GFF\n";
        read_fasta(\%SEQUENCE, $fasta_file);

        for my $trans (keys %SEQUENCE) {
            $trans =~ /^(.+$)(\s.*$|$)/;
            my $new_id = $1;
            if (grep { $new_id eq $_ } @ID_list) {
                my $for_print = ">$new_id\n$SEQUENCE{$trans}";
                push @for_print, $for_print;
            }
        }

        raid::MyFileIO::print_out(\@for_print, $output);
        my $seq_num = @for_print;
        print STDERR "==> Extract $seq_num sequences from GFF\n";
    }
}
else {
    print STDERR "==> No multi-transcripts.\n";

    # directly accept from fasta file
    read_fasta(\%SEQUENCE, $fasta_file);

    # write output
    for my $trans (keys %SEQUENCE) {
        $trans =~ /^(.+?)(\s.*$|$)/;
        my $new_id = $1;
        my $for_print = ">$new_id\n$SEQUENCE{$trans}";
        push @for_print, $for_print;
    }

    raid::MyFileIO::print_out(\@for_print, $output);
    my $seq_num = @for_print;
    print STDERR "==> All $seq_num sequences are still in.\n";
}


#print "Gene\tLongest_trans\n";

#for my $gene_id (keys %LONGEST) {
#    print "$gene_id\t$LONGEST{$gene_id}\n";
#}

#----------------------------------------------------------#
# sub
#----------------------------------------------------------#

# Usage : read_fasta(\%SEQUENCE, $fasta_file);
# return nothing but already read fasta into a hash

sub read_fasta{
    my ($hash, $in) = @_;

    my $seqIOobj = Bio::SeqIO -> new(
        -file       =>  "$in",
        '-format'   =>  'Fasta'
    );

    while((my $seqobj = $seqIOobj -> next_seq())) {
        my $id = $seqobj -> id();
        my $seq = $seqobj -> seq();
        $hash -> {$id} = $seq;
    }
}
