#!/usr/bin/env perl
#
# Convert.pm -- Multiple converting tools
#
# Author: Yuqian Jiang
# Created: 2023-06-08
# Version: 1.0.0
#
# Change logs:
# Version 1.0.0 2023-06-08: Initial version. Add function codon_translate
# Version 1.0.1 2023-06-15: conda_translate could translate sequences now

=head1 NAME

raid::Convert - Converting some results into the specific types

=head1 SYNOPSIS

    use raid::Convert qw();

=cut

package raid::Convert;
use strict;
use warnings;
use autodie;

=head1 METHODS

=head2 ECD_convert

=cut
sub ECD_convert {
    my @
}

=head2 codon_translate

      About : Tranlating CDS to proteins
      Usage : my $pep_seq = codon_translate($dna_seq)
       Args : codon translator
    Returns : pep sequences

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
    }
    my $out_pep = join '', @pep_seq;
    return $out_pep;
}


1;

=head1 VERSION

1.0.1

=head1 AUTHOR

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

Special thanks to Qiang Wang, wang-q@outlook.com
Special thanks to Nowind, noitulove9891@gmail.com

=cut
