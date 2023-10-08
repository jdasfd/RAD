#!/usr/bin/perl -w
#
# RAID.pl -- RLK Auto-IDentifier
#
# Author: Yuqian Jiang
# Created: 2023-08-10
# Version: 1.2.4
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
use raid::MyFileIO;
use raid::DomUtil;
use raid::OptSeq;


#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

RAID.pl -- RLK Auto-IDentifier

=head1 SYNOPSIS

    RAID.pl (version 1.2.4)
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
my $i = 0;
my $gene_name;
my $tee_add = IO::Tee -> new ( ">> $outdir/log.txt", \*STDERR );
my (%DOMAIN_info, %SEQUENCE, %TM_info, %RLK);
my @PKD = ("Pkinase", "PK_Tyr_Ser-Thr", "Pkinase_fungal", "Pkinase_C");
my (@query_with_pkd, @write_hmm_tsv, @final_domain_tsv, @rlk_out_tsv, @other_rlk_out_tsv);

# sequence read to hash
# start information record to log
my @all_seq_id = raid::MyFileIO::read_fasta(\%SEQUENCE, $input);
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
my @hmmresult = raid::MyFileIO::read_hmm_txt($pfam_txt);

# input hmm result to a hash
for ( @hmmresult ) {
    my @col = split/\t/, $_;
    if ( !defined $gene_name ) {
        $gene_name = $col[0];
        my $info = "$col[1],$col[2],$col[3],$col[4]";
        $DOMAIN_info_ref -> {$col[0]} -> [$i] = $info;
        $i++;
    }
    elsif ( $gene_name eq $col[0] ) {
        my $info = "$col[1],$col[2],$col[3],$col[4]";
        $DOMAIN_info_ref -> {$col[0]} -> [$i] = $info;
        $i++;
    }
    elsif ( $gene_name ne $col[0] ) {
        $gene_name = $col[0];
        $i = 0;
        my $info = "$col[1],$col[2],$col[3],$col[4]";
        $DOMAIN_info_ref -> {$col[0]} -> [$i] = $info;
        $i++;
    }
}


#----------------------------------------------------------#
# Kinase domain extraction
#----------------------------------------------------------#

# sort hash ascendingly according to e-value
# domain_sort function could be seen in the perl module locally
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

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
raid::MyFileIO::print_out(\@write_hmm_tsv, $pfam_tsv);

# combined with pkd results.
@query_with_pkd = uniq (@query_with_pkd);
my $kd_num = @query_with_pkd;
if ( $kd_num != 0 ) {
    print $tee_add "Detected totally $kd_num proteins with Pkinase domains.\n";
    print $tee_add "\n";
}
else {
    print $tee_add "None Pkinase domains detected. Aborted!\n";
    die;
}

# Proteins with KD saved to list and fa files
raid::MyFileIO::print_out(\@query_with_pkd, $pro_with_kd);

my @KD_fa = raid::OptSeq::seq_some(\%SEQUENCE, \@query_with_pkd);
raid::MyFileIO::print_out(\@KD_fa, $kd_seq);


#----------------------------------------------------------#
# TMbed
#----------------------------------------------------------#

# tmbed
if ( !path($tmbed_pred) -> is_file || defined $forcetmd ) {
    print $tee_add "==> TMbed start.\n";
    my $tmbedresult = tmbed($kd_seq, $batchsize, $tmbed_pred);
    if ($tmbedresult != 0) {
        print $tee_add "Error: tmbed interrupted.\n";
        die;
    }
}
else {
    print $tee_add "==> TMbed results already exists. Moving to the next step\n";
}

# dealing with the tmbed result
raid::MyFileIO::read_pred(\%TM_info, $tmbed_pred);
my $pred_num = keys %TM_info;
print $tee_add "TMbed result contains $pred_num proteins.\n";
if ( $kd_num == $pred_num ) {
    print $tee_add "All proteins with pkinase were searched via TMbed.\n";
}
else {
    print $tee_add "Something went wrong when running tmbed, please check tmbed sequences.\n";
    die;
}

raid::MyFileIO::extract_pred_info(\%DOMAIN_info, \%TM_info);
# sort all domains with SP and TMD inside
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "2", ",");

#print Dumper \%DOMAIN_info;
print $tee_add "\n";


#----------------------------------------------------------#
# Domain manipulating
#----------------------------------------------------------#

print $tee_add "==> Start to deal with domains.\n";

# domain filter and sort according to start pos
raid::DomUtil::domain_filter(\%DOMAIN_info);
raid::DomUtil::domain_sort(\%DOMAIN_info, "4", "3", ",");

for ( @query_with_pkd ) {
    my $id = $_;
    my $array_ref = $DOMAIN_info_ref -> {$id};

    for ( @{$array_ref} ) {
        my @output = split /,/, $_;
        if ( grep { $_ eq $output[0] } @PKD ) {
            $output[0] = "Kinase";
        }
        my $outline = "$id\t$output[0]\t$output[2]\t$output[3]";
        push @final_domain_tsv, $outline;
    }
}

# write to final tsv (only contained all pkinase pros)
print $tee_add "\n";
print $tee_add "==> Final domain result output.\n";
raid::MyFileIO::print_out(\@final_domain_tsv, $domain_final);


#----------------------------------------------------------#
# RLK scanning
#----------------------------------------------------------#

my $headline = "Name\tType\tECD\tKD_count";
push @rlk_out_tsv, $headline;
push @other_rlk_out_tsv, $headline;
print $tee_add "\n";
print $tee_add "==> Scanning RLKs.\n";

# avoid confusing, reinduced the final domain tsv
my $final_in = raid::MyFileIO::getInputFilehandle($domain_final);
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
            $ECD_all =~ s/TMD_[i|o]2[o|i]//g;
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
            my $outline = "$keys\tRLK_reverse\tNone\t$KD_count";
            push @other_rlk_out_tsv, $outline;
        }
        else {
            my $outline = "$keys\tRLK_WE\tNone\t$KD_count";
            push @rlk_out_tsv, $outline;
        }
    }
    else {
        my $outline = "$keys\tOthers\tNone\t$KD_count";
        push @other_rlk_out_tsv, $outline;
    }
}

raid::MyFileIO::print_out(\@rlk_out_tsv, $rlk_output);
raid::MyFileIO::print_out(\@other_rlk_out_tsv, $other_rlk_output);

my $rlk_count_num = @rlk_out_tsv;
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

1.2.4

=head1 AUTHORS

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
