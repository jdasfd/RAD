#!/usr/bin/env perl
#
# OptSeq.pm -- Multiple sequences operation
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.2.4
#
# Change logs:
# Version 1.0.0 2023-06-08: Initial version. Add function codon_translate.
# Version 1.0.1 2023-06-15: Codon_translate could translate sequences now.
# Version 1.1.0 2023-08-13: Rename .pm to OptSeq.
# Version 1.1.1 2023-08-14: Add function seq_some, seq_replace.
# Version 1.2.0 2023-08-14: Add function seq_trunc, but not finished yet.
# Version 1.2.1 2023-08-15: Bug fixes for fasta output without > as id start.
# Version 1.2.2 2024-01-08: Complete seq_trunc.
# Version 1.2.3 2024-06-13: Add X to represent the abnormal cds.
# Version 1.2.4 2025-04-22: Change package name to rad

=head1 NAME

rad::OptSeq - Converting some results into the specific types

=head1 SYNOPSIS

    use rad::OptSeq qw();

=cut

package rad::OptSeq;
use strict;
use warnings;
use autodie;

=head1 METHODS

=head2 codon_translate

      About : Tranlating CDS to proteins.
      Usage : my $pep_seq = codon_translate($dna_seq)
       Args : Codon translator.
    Returns : Pep sequences.

=cut
sub codon_translate {
    my ($in_cds) = shift;

    my (%codon) = (
        'TCA' => 'S',   # Serine
        'TCC' => 'S',   # Serine
        'TCG' => 'S',   # Serine
        'TCT' => 'S',   # Serine
        'TTC' => 'F',   # Phenylalanine
        'TTT' => 'F',   # Phenylalanine
        'TTA' => 'L',   # Leucine
        'TTG' => 'L',   # Leucine
        'TAC' => 'Y',   # Tyrosine
        'TAT' => 'Y',   # Tyrosine
        'TAA' => '*',   # Stop
        'TAG' => '*',   # Stop
        'TGC' => 'C',   # Cysteine
        'TGT' => 'C',   # Cysteine
        'TGA' => '*',   # Stop
        'TGG' => 'W',   # Tryptophan
        'CTA' => 'L',   # Leucine
        'CTC' => 'L',   # Leucine
        'CTG' => 'L',   # Leucine
        'CTT' => 'L',   # Leucine
        'CCA' => 'P',   # Proline
        'CCC' => 'P',   # Proline
        'CCG' => 'P',   # Proline
        'CCT' => 'P',   # Proline
        'CAC' => 'H',   # Histidine
        'CAT' => 'H',   # Histidine
        'CAA' => 'Q',   # Glutamine
        'CAG' => 'Q',   # Glutamine
        'CGA' => 'R',   # Arginine
        'CGC' => 'R',   # Arginine
        'CGG' => 'R',   # Arginine
        'CGT' => 'R',   # Arginine
        'ATA' => 'I',   # Isoleucine
        'ATC' => 'I',   # Isoleucine
        'ATT' => 'I',   # Isoleucine
        'ATG' => 'M',   # Methionine (start codon)
        'ACA' => 'T',   # Threonine
        'ACC' => 'T',   # Threonine
        'ACG' => 'T',   # Threonine
        'ACT' => 'T',   # Threonine
        'AAC' => 'N',   # Asparagine
        'AAT' => 'N',   # Asparagine
        'AAA' => 'K',   # Lysine
        'AAG' => 'K',   # Lysine
        'AGC' => 'S',   # Serine
        'AGT' => 'S',   # Serine
        'AGA' => 'R',   # Arginine
        'AGG' => 'R',   # Arginine
        'GTA' => 'V',   # Valine
        'GTC' => 'V',   # Valine
        'GTG' => 'V',   # Valine
        'GTT' => 'V',   # Valine
        'GCA' => 'A',   # Alanine
        'GCC' => 'A',   # Alanine
        'GCG' => 'A',   # Alanine
        'GCT' => 'A',   # Alanine
        'GAC' => 'D',   # Aspartic Acid
        'GAT' => 'D',   # Aspartic Acid
        'GAA' => 'E',   # Glutamic Acid
        'GAG' => 'E',   # Glutamic Acid
        'GGA' => 'G',   # Glycine
        'GGC' => 'G',   # Glycine
        'GGG' => 'G',   # Glycine
        'GGT' => 'G',   # Glycine
    );

    my $in_len = length($in_cds);
    unless ( $in_len%3 == 0 ) {
        print STDERR "Warning: input sequence is not a interger multiple of 3."
    }

    $in_cds = uc($in_cds);
    my @nt_seq = unpack("(A3)*", $in_cds);
    my @pep_seq = ();
    for ( my $i = 0; $i < @nt_seq; $i++ ) {
        if ($codon{$nt_seq[$i]}) {
            push @pep_seq, $codon{$nt_seq[$i]};
        }
        else {
            push @pep_seq, "X";
        }
    }
    my $out_pep = join '', @pep_seq;
    return $out_pep;
}

=head2 seq_some

      About : Extracting sequences from sequence hash.
      Usage : my @for_print = seq_some(\%SEQUENCE, \@some_id);
       Args : Hash with all sequences;
              id list for extracting sequences.
    Returns : Printing format of sequences.

=cut
sub seq_some {
    my ($seq_hash_ref, $id_array_ref) = @_;
    my @for_print;

    for ( @{$id_array_ref} ) {
        my $seq = $seq_hash_ref -> {$_};
        push @for_print, ">$_";
        push @for_print, $seq;
    }

    return @for_print;
}

=head2 seq_replace

      About : Replacing sequence names.
      Usage : my @for_print = seq_some(\%SEQUENCE, \@replace_id);
       Args : Hash with all sequences;
              replacing array with element saved "old_name,new_name".
    Returns : Printing format of sequences.

=cut
sub seq_replace {
    my ($seq_hash_ref, $rep_array_ref) = @_;
    my @for_print;

    for ( @{$rep_array_ref} ) {
        my ($old_id, $new_id) = split/,/, $_;
        my $seq = $seq_hash_ref -> {$old_id};
        push @for_print, $new_id;
        push @for_print, $seq;
    }

    return @for_print;
}

=head2 seq_trunc

      About : Truncating sequences according to domain info.
      Usage : my @for_print = seq_trunc(\%SEQUENCE, \%SEQ_INFO, \@id);
       Args : Hash with all sequences;
              hash with truncated domain info (domain,start,stop without e-value info);
              array saved all id.
    Returns : Printing format of sequences.

=cut
sub seq_trunc {
    my ($seq_ref, $seqinfo_ref, $id_ref) = @_;
    my @for_print;

    for ( @{$id_ref} ) {
        my $id = $_;
        my $seqinfo = $seqinfo_ref -> {$id};
        for ( @{$seqinfo} ) {
            my ($domain, $start, $stop) = split/,/, $_;
            my $seq_name = "$id"."_"."$domain"."_"."$start"."_"."$stop";
            my $sequence = $seq_ref -> {$id};
            my $seq_start = $start - 1;
            my $seq_len = $stop - $start + 1;
            my $truc_seq = substr($sequence, $seq_start, $seq_len);
            push @for_print, ">$seq_name\n$truc_seq";
        }
    }

    return @for_print;
}

1;

=head1 VERSION

1.2.4

=head1 AUTHOR

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

Special thanks to Qiang Wang, wang-q@outlook.com
Special thanks to Nowind, noitulove9891@gmail.com

=cut
