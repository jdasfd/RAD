#!/usr/bin/perl -w
#
# RLPscanner.pl -- RLP Auto-Identifier
#
# Author: Yuqian Jiang
# Created: 2024-04-06

# Change logs:
# 2024-04-06: The initial version. Realize automatically RLP scanning.

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use List::Util qw(uniq);
use Data::Dumper;
use Math::BigFloat;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use raid::MyFileIO;
use raid::DomUtil;
use raid::OptSeq;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

RLPscanner.pl -- RLP Auto-IDentifier

=head1 SYNOPSIS

    RLPscanner.pl
        RLP Automatical IDentifier searching RLPs among protein.fa files.

    Usage:

    Required:
        -d,--domain        STR       input domain scanning results, required
        -t,--tmbed         STR       transmembrane domain scanning results, required

    Options:
        -o,--outdir        STR       output dir
        -h,--help                    help information

=cut

GetOptions(
    "d|domain=s"        => \(my $domain),
    "t|tmbed=s"         => \(my $tmbed),
    "o|outdir=s"        => \(my $outdir),
    "h|help"            => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $domain ) {
    print STDERR "Error: cannot find input files.\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $tmbed ) {
    print STDERR "Error: cannot find hmm directory. Please check your path\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $outdir ) {
    my $current = Path::Tiny -> cwd;
    $outdir = $current."/result";
}
else {
    if ( $outdir =~ /\/$/ ) {
        $outdir =~ s/\/$//;
    }
    path ($outdir) -> mkdir;
}

#----------------------------------------------------------#
# Prep
#----------------------------------------------------------#

# value usage
my $defaulte = Math::BigFloat -> new('1e-3');
my (%DOMAIN_info, %TM_info, %RLP_like);
my (@tmdlist, @final_domain_tsv, @rlp_out_tsv);

# all output path
my $domain_final = $outdir."/Pro.final.domain.tsv";
my $rlp_output = $outdir."/RLP.tsv";

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#

# input the scan result via hmmscan
my $DOMAIN_info_ref = \%DOMAIN_info;
raid::MyFileIO::read_hmm_txt($DOMAIN_info_ref, $domain);
# domain_sort function could be seen in the perl module locally
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

# dealing with the tmbed result
raid::MyFileIO::read_pred(\%TM_info, $tmbed);

raid::MyFileIO::extract_pred_info(\%DOMAIN_info, \%TM_info);
# sort all domains with SP and TMD inside
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

# domain filter and sort according to start pos
raid::DomUtil::domain_filter(\%DOMAIN_info);
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "3", ",");

# get all pros with tmd
for my $gene ( keys %DOMAIN_info ) {
    my $info = $DOMAIN_info_ref -> {$gene};
    foreach ( @{$info} ) {
        my @array = split/,/, $_;
        if ( $array[0] =~ /TMD/ ) {
            push @tmdlist, $gene;
        }
    }
}

for ( @tmdlist ) {
    my $id = $_;
    my $domain_ref = $DOMAIN_info_ref -> {$id};

    for ( @{$domain_ref} ) {
        my @output = split /,/, $_;
        my $dom_array_e = Math::BigFloat -> new ($output[1]);
        if ( $dom_array_e <= $defaulte ) {
            my $outline = "$id\t$output[0]\t$output[2]\t$output[3]";
            push @final_domain_tsv, $outline;
        }
    }
}

# write to final tsv with TMD and SP
raid::MyFileIO::print_out(\@final_domain_tsv, $domain_final);

my $headline = "Name\tType\tECD";
push @rlp_out_tsv, $headline;

# reinduced the final domain tsv
my $final_in = raid::MyFileIO::getInputFilehandle($domain_final);
while ( <$final_in> ) {
    chomp;
    my @array = split/\t/, $_;
    push @{$RLP_like{$array[0]}}, $array[1];
}

my $domain_list;
# define and check domain structure
for my $keys (keys %RLP_like) {
    $domain_list = join("#", @{$RLP_like{$keys}});
    if ( $domain_list =~ /TMD_/ ) {
        if ($domain_list =~ /^(.+?)TMD_[i|o]2[o|i]$/) {
            my $ECD_all = $1;
            $ECD_all =~ s/^Sig_Pep#//g;
            $ECD_all =~ s/TMD_[i|o]2[o|i]#//g;
            $ECD_all =~ s/#$//;
            if ($ECD_all ne "") {
                my $outline = "$keys\tRLP\t$ECD_all";
                push @rlp_out_tsv, $outline;
            }
            else {
                my $outline = "$keys\tRLPUN\tUnknown";
                push @rlp_out_tsv, $outline;
            }
        }
    }
}

# write into rlp tsv
raid::MyFileIO::print_out(\@rlp_out_tsv, $rlp_output);
