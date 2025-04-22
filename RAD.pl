#!/usr/bin/perl -w
#
# RAD.pl -- RLK Auto-IDentifier
#
# Author: Yuqian Jiang
# Created: 2023-08-10
# Version: 1.2.8
#
# Change logs:
# Version 1.0.0 2023-08-10: The initial version. Realize automatically RLK scanning from a protein fasta file.
# Version 1.1.0 2023-08-15: All tested in my own environment and run successfully. Remove repeated modules.
# Version 1.1.1 2023-08-17: Resolve log cannot redirect to the specified output.
#                           Remove File::Basename for not using it.
#                           Create output dir.
# Version 1.2.0 2023-09-28: Judge different kinases - will judge the first kinase after the first TMD.
#                           Log will add system time and will not replaced old log.
#                           RLK judgement standard change, modified code of RLK scanning part.
# Version 1.2.1 2023-09-28: RLK number counting added to log. ECD output modified (remove the last -)
# Version 1.2.2 2023-10-06: More TMD RLKs would not be removed.
# Version 1.2.3 2023-10-07: Bug fixes: TMD_i2o would be treated as ECD in some cases.
# Version 1.2.4 2023-10-08: All ECD were unified to None for better viewing. Add headline to RLK.others.tsv.
#                           Bug fixes: RLK numbers count without headline.
# Version 1.2.5 2023-11-13: Bug fixes: Some RLKs will get only # as ECD (which means they are RLK_WE)
# Version 1.2.6 2024-01-08: Bug fixes: removing those fake RLKs and probably fake domains.
# Version 1.2.7 2024-01-23: Bug fixes: RLKs finding process. Will count something not RLKs.
# Version 1.2.8 2025-04-22: Change name to RAD.pl

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use List::Util qw(uniq);
use IO::Tee;
use Data::Dumper;
use Math::BigFloat;
use FindBin qw/$Bin/;
use lib "$FindBin::Bin/lib/";
use rad::MyFileIO;
use rad::DomUtil;
use rad::OptSeq;


#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

RAD.pl -- RLK Auto-IDentifier

=head1 SYNOPSIS

    RAD.pl (version 1.2.8)
        RLK Automatical IDentifier searching RLKs among protein.fa files.

    Usage:

    Required:
        -i,--in          STR       input protein.fa, required
        --hmm            STR       hmm database for hmmscan, required

    Options:
        -o,--outdir      STR       out dir, default is creating a dir named result under the current path
        -e,--evalue      FLO       E-value for hmmscan, default: 0.1 (see -E and --domE in hmmscan -h)
        --thread         INT       threads for hmmscan, default: 1
        --forcescan                force hmmscan run again, default: skip hmmscan if the result exists
        --batchsize      INT       batch size for tmbed usage, default: 1
        --forcetmd                 force tmbed run again, default: skip tmbed if the result exists
        -h,--help                  help information

=cut

GetOptions(
    "i|in=s"            => \(my $input),
    "hmm=s"             => \(my $hmmpath),
    "o|outdir=s"        => \(my $outdir),
    "e|evalue=f"        => \(my $e_value = '0.1'),
    "thread=i"          => \(my $thread = '1'),
    "forcescan"         => \(my $forcescan),
    "batchsize=i"       => \(my $batchsize = '1'),
    "forcetmd"          => \(my $forcetmd),
    "h|help"            => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: cannot find input files.\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $hmmpath ) {
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

unless ( path("$outdir/log.txt") -> is_file ) {
    path("$outdir/log.txt") -> touch;
}


#----------------------------------------------------------#
# Prep
#----------------------------------------------------------#

# value usage
my $defaulte = Math::BigFloat -> new('1e-10');
my $defaulte_others = Math::BigFloat -> new('1e-3');
my $tee_add = IO::Tee -> new ( ">> $outdir/log.txt", \*STDERR );
my (%DOMAIN_info, %SEQUENCE, %TM_info, %RLK);
my @PKD = ("Pkinase", "PK_Tyr_Ser-Thr", "Pkinase_fungal", "Pkinase_C");
my (@query_with_pkd, @write_hmm_tsv, @final_domain_tsv, @rlk_out_tsv, @other_rlk_out_tsv);

# sequence read to hash
# start information record to log
my @all_seq_id = rad::MyFileIO::read_fasta(\%SEQUENCE, $input);
my $all_seq_num = @all_seq_id;
print $tee_add "Start to scan $all_seq_num proteins inside the fasta file!\n";
print $tee_add "Input: $input.\n";
my $time = localtime;
print $tee_add "Start at $time.\n";
print $tee_add "\n";

# all output path
my $pfam_txt = $outdir."/pfam.txt";
my $pfam_tsv = $outdir."/pfam.tsv";
my $pro_with_kd = $outdir."/Pro.KD.lst";
my $kd_seq = $outdir."/Pro.KD.fa";
my $tmbed_pred = $outdir."/Pro.KD.pred";
my $domain_final = $outdir."/Pro.final.domain.tsv";
my $rlk_output = $outdir."/RLK.tsv";
my $other_rlk_output = $outdir."/RLK.others.tsv";


#----------------------------------------------------------#
# HMMSCAN
#----------------------------------------------------------#

# starting with hmmscan
if ( !path($pfam_txt) -> is_file || defined $forcescan ) {
    print $tee_add "==> HMMSCAN start.\n";
    my $hmmresult = hmmscan($input, $thread, $e_value, $outdir, $hmmpath);
    if ($hmmresult != 0) {
        print $tee_add "Error: hmmscan interrupted.\n";
        die;
    }
}
else {
    print $tee_add "==> HMMSCAN results already exists. Moving to the next step\n";
}

# read the scan result via hmmscan
my $DOMAIN_info_ref = \%DOMAIN_info;
rad::MyFileIO::read_hmm_txt($DOMAIN_info_ref, $pfam_txt);


#----------------------------------------------------------#
# Kinase domain extraction
#----------------------------------------------------------#

# sort hash ascendingly according to e-value
# domain_sort function could be seen in the perl module locally
rad::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

for my $gene ( keys %DOMAIN_info ) {
    my $info = $DOMAIN_info_ref -> {$gene};
    foreach ( @{$info} ) {
        my @array = split/,/, $_;
        my $line = "$gene\t$array[0]\t$array[1]\t$array[2]\t$array[3]";
        push @write_hmm_tsv, $line;
        my $array_e = Math::BigFloat -> new ($array[1]);
        if ( grep { $_ eq $array[0] } @PKD ) {
            if ( $array_e <= $defaulte ) {
                push @query_with_pkd, $gene;
            }
            else {
                next;
            }
        }
    }
}

# output into tsv according to the hmm txt result
rad::MyFileIO::print_out(\@write_hmm_tsv, $pfam_tsv);

# combined with pkd results.
@query_with_pkd = uniq (@query_with_pkd);
my $kd_num = @query_with_pkd;
if ( $kd_num != 0 ) {
    print $tee_add "Detected totally $kd_num proteins with true Pkinase domains.\n";
    print $tee_add "\n";
}
else {
    print $tee_add "None Pkinase domains detected. Aborted!\n";
    die;
}

# Proteins with KD saved to list and fa files
rad::MyFileIO::print_out(\@query_with_pkd, $pro_with_kd);

my @KD_fa = rad::OptSeq::seq_some(\%SEQUENCE, \@query_with_pkd);
rad::MyFileIO::print_out(\@KD_fa, $kd_seq);


#----------------------------------------------------------#
# TMbed
#----------------------------------------------------------#

# tmbed
if ( !path($tmbed_pred) -> is_file || defined $forcetmd ) {
    print $tee_add "==> TMbed start.\n";
    my $tmbedresult = tmbed($kd_seq, $batchsize, $tmbed_pred);
    if ($tmbedresult != 0) {
        print $tee_add "Error: tmbed interrupted.\n";
        print $tee_add "\n";
        die;
    }
}
else {
    print $tee_add "==> TMbed results already exists. Moving to the next step\n";
}

# dealing with the tmbed result
rad::MyFileIO::read_pred(\%TM_info, $tmbed_pred);
my $pred_num = keys %TM_info;
print $tee_add "TMbed result contains $pred_num proteins.\n";
if ( $kd_num == $pred_num ) {
    print $tee_add "All proteins with pkinase were searched via TMbed.\n";
}
else {
    print $tee_add "Something went wrong when running tmbed, please check tmbed sequences.\n";
    print "\n";
    die;
}

rad::MyFileIO::extract_pred_info(\%DOMAIN_info, \%TM_info);
# sort all domains with SP and TMD inside
rad::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

#print Dumper \%DOMAIN_info;
print $tee_add "\n";


#----------------------------------------------------------#
# Domain manipulating
#----------------------------------------------------------#

print $tee_add "==> Start to deal with domains.\n";

# domain filter and sort according to start pos
rad::DomUtil::domain_filter(\%DOMAIN_info);
rad::DomUtil::domain_sort(\%DOMAIN_info, "4", "3", ",");

for ( @query_with_pkd ) {
    my $id = $_;
    my $domain_ref = $DOMAIN_info_ref -> {$id};

    for ( @{$domain_ref} ) {
        my @output = split /,/, $_;
        my $dom_array_e = Math::BigFloat -> new ($output[1]);
        if ( grep { $_ eq $output[0] } @PKD ) {
            if ( $dom_array_e <= $defaulte ) {
                $output[0] = "Kinase";
            }
            else {
                $output[0] = "Fake_kinase";
            }
            my $outline = "$id\t$output[0]\t$output[2]\t$output[3]";
            push @final_domain_tsv, $outline;
        }
        else {
            if ( $dom_array_e <= $defaulte_others ) {
                my $outline = "$id\t$output[0]\t$output[2]\t$output[3]";
                push @final_domain_tsv, $outline;
            }
        }
    }
}

# write to final tsv (only contained all pkinase pros)
print $tee_add "\n";
print $tee_add "==> Final domain result output.\n";
rad::MyFileIO::print_out(\@final_domain_tsv, $domain_final);


#----------------------------------------------------------#
# RLK scanning
#----------------------------------------------------------#

my $headline = "Name\tType\tECD\tKD_count";
push @rlk_out_tsv, $headline;
push @other_rlk_out_tsv, $headline;
print $tee_add "\n";
print $tee_add "==> Scanning RLKs.\n";

# avoid confusing, reinduced the final domain tsv
my $final_in = rad::MyFileIO::getInputFilehandle($domain_final);
while ( <$final_in> ) {
    chomp;
    my @array = split/\t/, $_;
    push @{$RLK{$array[0]}}, $array[1];
}

my $domain_list;
# define and check domain structure
for my $keys (keys %RLK) {
    $domain_list = join("#", @{$RLK{$keys}});
    my $KD_count = &COUNT_SUB_STR("Kinase", $domain_list);
    if ( $domain_list =~ /TMD_/ ) {
        if ($domain_list =~ /^(.+?)TMD_o2i.+?Kinase.*$/) {
            my $ECD_all = $1;
            $ECD_all =~ s/^Sig_Pep#//g;
            $ECD_all =~ s/TMD_[i|o]2[o|i]#//g;
            $ECD_all =~ s/#$//;
            if ($ECD_all ne "") {
                my $outline = "$keys\tRLK\t$ECD_all\t$KD_count";
                push @rlk_out_tsv, $outline;
            }
            else {
                my $outline = "$keys\tRLK_WE\tNone\t$KD_count";
                push @rlk_out_tsv, $outline;
            }
        }
        elsif ($domain_list =~ /TMD_i2o/) {
            my $outline = "$keys\tRLK_RvrTMD\tNone\t$KD_count";
            push @other_rlk_out_tsv, $outline;
        }
        else {
            if ( $domain_list =~ /^TMD_o2i.+?Kinase.*/ ) {
                my $outline = "$keys\tRLK_WE\tNone\t$KD_count";
                push @rlk_out_tsv, $outline;
            }
            else {
                my $outline = "$keys\tOthers\tNone\t$KD_count";
                push @other_rlk_out_tsv, $outline;
            }
        }
    }
    else {
        my $outline = "$keys\tOthers\tNone\t$KD_count";
        push @other_rlk_out_tsv, $outline;
    }
}

rad::MyFileIO::print_out(\@rlk_out_tsv, $rlk_output);
rad::MyFileIO::print_out(\@other_rlk_out_tsv, $other_rlk_output);


#----------------------------------------------------------#
# Finish
#----------------------------------------------------------#

my $rlk_count_num = @rlk_out_tsv;
$rlk_count_num = $rlk_count_num - 1;
print $tee_add "$rlk_count_num RLKs scanned.\n";
my $end_time = localtime;
print $tee_add "==> Finished at $end_time!\n";
print $tee_add "\n";


#----------------------------------------------------------#
# Sub program
#----------------------------------------------------------#

# Usage: my $result = hmmscan($input, $thread, $e_value, $out, $hmmpath);
sub hmmscan {
    my ($INPATH, $THREAD, $E_VALUE, $OUT, $HMMPATH) = @_;

    my $result = system "hmmscan --cpu $THREAD -E $E_VALUE --domE $E_VALUE --noali --notextw -o $OUT/pfam.txt $HMMPATH $INPATH > /dev/null";
    return $result;
}

sub tmbed {
    my ($INPATH, $BATCH, $OUT) = @_;

    my $result = system "tmbed predict --batch-size $BATCH --use-gpu -f $INPATH -p $OUT";
    return $result;
}

sub COUNT_SUB_STR {
    my $SUB_STR = $_[0];
    my $STR = $_[1];
    # count sub str
    my $tmp = 0;
    my $count = 0;
    if ( $tmp = () = ($STR) =~ /$SUB_STR/g ) {
        $count += $tmp;
    }
    return $count;
}

=head1 VERSION

1.2.8

=head1 AUTHORS

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
