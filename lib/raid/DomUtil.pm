#!/usr/bin/env perl
#
# DomUtil.pm -- Domain manipulating tools
#
# Author: Yuqian Jiang
# Created: 2023-07-20
# Version: 1.2.0
#
# Change logs:
# Version 1.0.0 2023-07-20: Initial version. Add function domain_unified but need updating.
# Version 1.1.0 2023-08-13: Add function domain_sort.
# Version 1.2.0 2023-08-15: Add function domain_filter.

=head1 NAME

raid::DomUtil - Local perl module for operating domains

=head1 SYNOPSIS

    use raid::DomUtil;

=cut

package raid::DomUtil;

use strict;
use warnings;
use AlignDB::IntSpan;
use Sort::Array qw (Sort_Table);

=head1 METHODS

=head2 domain_unified

      About : Converting different domains into unified names
      Usage : my $unified_name = domain_unified($ECD_name);
       Args : codon translator
    Returns : domain name unified

=cut
sub ECD_convert {
    my $ECD = @_;
}

=head2 domain_sort

      About : Sort domains
      Usage : domain_sort(\%DOMAIN, $all_col, $sorted_col, $sep);
       Args : hash saved all domain info, sort inside value
    Returns : Nothing

=cut
sub domain_sort {
    my ($seq_hash, $allc, $sortc, $sep) = @_ ;

    for my $gene ( keys %{$seq_hash} ) {
        my $dom_info = $seq_hash -> {$gene};
        my @domain_for_sort = @{$dom_info};
        my @domain_sorted = Sort_Table(
            cols      => $allc,
            field     => $sortc,
            sorting   => 'ascending',
            structure => 'csv',
            separator => $sep,
            data      => \@domain_for_sort,
        );
        @{$dom_info} = @domain_sorted;
        $seq_hash -> {$gene} = $dom_info;
    }
}

=head2 domain_filter

      About : Filter all sorted domains according to their positions.
      Usage : domain_filter();
       Args : Hash saved all domain info, sort with e-value from less to more.
    Returns : Nothing.

=cut
sub domain_filter {
    my ($domain_hash) = @_;

    for my $gene_id ( keys %{$domain_hash} ) {
        my $domain_array_ref = $domain_hash -> {$gene_id};
        my $domain_set = AlignDB::IntSpan -> new;
        my @new_domain_list;

        for ( @{$domain_array_ref} ) {
            my $record = $_;
            my ($start, $end) = (split /,/, $_) [2,3];
            if ( $domain_set -> is_empty ) {
                $domain_set -> add_range($start, $end);
                push @new_domain_list, $record;
            }
            else {
                my $test_set = AlignDB::IntSpan -> new;
                $test_set -> add_range($start, $end);
                my $result = $domain_set -> intersect($test_set);
                if ( $result -> is_empty ) {
                    $domain_set -> add_range($start, $end);
                    push @new_domain_list, $record;
                }
                else {
                    next;
                }
            }
        }

        @{$domain_array_ref} = @new_domain_list;
        $domain_hash -> {$gene_id} = $domain_array_ref;
    }
}

1;

=head1 VERSION

1.2.0

=head1 AUTHOR

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
