#!/usr/bin/perl

# Copyright 2010 Andreas R. Gruber
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details (http://www.gnu.org/licenses).

use strict;
use warnings;

###############################################################################
# MAIN
###############################################################################

my $motif = $ARGV[0];
my $file  = $ARGV[1];

exit 0 if ( !-e $file );

open( FILE, $file );
my $tmp = '';

my $overall_counts   = 0;
my $theo_counts_mono = 0;
my $theo_counts_di   = 0;

while ( my $line = <FILE> ) {
  chomp $line;
  if ( $line =~ m/^>/ ) {
    if ( $tmp ne '' ) {
      my $nuc_counts = _count_nucleotides($tmp);
      my $counts = _find_motif( $tmp, $motif );
      my ( $mono, $di ) = _theoretical_occurences( $nuc_counts, $motif );
      $overall_counts   += $counts;
      $theo_counts_mono += $mono;
      $theo_counts_di   += $di;
    }
    $tmp = '';
    next;
  }
  $tmp .= $line;
}

my $nuc_counts = _count_nucleotides($tmp);
my $counts = _find_motif( $tmp, $motif );
my ( $mono, $di ) = _theoretical_occurences( $nuc_counts, $motif );
$overall_counts   += $counts;
$theo_counts_mono += $mono;
$theo_counts_di   += $di;

print "                     MOTIF: $motif\n";
print "                      FILE: $file\n";
print "                    COUNTS: $overall_counts\n";
print "THEORETICAL MONONUCLEOTIDE: $theo_counts_mono\n";
print "  THEORETICAL DINUCLEOTIDE: $theo_counts_di\n";

###############################################################################
# SUBS
###############################################################################

sub _theoretical_occurences {
  my $nuc_counts = $_[0];
  my $motif      = $_[1];

  my $p_mono = 1;
  my $p_di   = 1;

  for my $i ( 0 .. length($motif) - 2 ) {
    my $nuc   = substr( $motif, $i, 1 );
    my $dinuc = substr( $motif, $i, 2 );

    if ( defined $nuc_counts->{$nuc} ) {
      $p_mono = $p_mono * ( $nuc_counts->{$nuc} / $nuc_counts->{'LENGTH'} );
    } else {
      $p_mono = 0;
    }
    if ( defined $nuc_counts->{$dinuc} ) {
      $p_di = $p_mono * ( $nuc_counts->{$dinuc} / $nuc_counts->{'LENGTH'} );
    } else {
      $p_di = 0;
    }
  }
  my $lastnuc = substr( $motif, length($motif) - 1, 1 );
  if ( defined $nuc_counts->{$lastnuc} ) {
    $p_mono = $p_mono * ( $nuc_counts->{$lastnuc} / $nuc_counts->{'LENGTH'} );
  } else {
    $p_mono = 0;
  }
  my $theo_counts_mono = sprintf( "%.5f", $nuc_counts->{'LENGTH'} * $p_mono );
  my $theo_counts_di   = sprintf( "%.5f", $nuc_counts->{'LENGTH'} * $p_di );

  return ( $theo_counts_mono, $theo_counts_di );
}

sub _find_motif {
  my $seq   = $_[0];
  my $motif = $_[1];

  $seq =~ s/\s//g;
  $seq = uc($seq);

  my $l      = length($motif);
  my $counts = 0;
  for my $i ( 0 .. length($seq) - $l ) {
    my $tmp = substr( $seq, $i, $l );
    $counts++ if ( $tmp eq $motif );
  }

  return $counts;
}

sub _count_nucleotides {
  my $seq = $_[0];

  $seq =~ s/\s//g;
  $seq = uc($seq);

  my $hash = {
    'A'      => 0,
    'C'      => 0,
    'G'      => 0,
    'T'      => 0,
    'AA'     => 0,
    'AC'     => 0,
    'AG'     => 0,
    'AT'     => 0,
    'CA'     => 0,
    'CC'     => 0,
    'CG'     => 0,
    'CT'     => 0,
    'GA'     => 0,
    'GC'     => 0,
    'GG'     => 0,
    'GT'     => 0,
    'TA'     => 0,
    'TC'     => 0,
    'TG'     => 0,
    'TT'     => 0,
    'LENGTH' => length($seq)
  };

  for my $i ( 0 .. $hash->{'LENGTH'} - 2 ) {
    my $nuc   = substr( $seq, $i, 1 );
    my $dinuc = substr( $seq, $i, 2 );
    $hash->{$nuc}++;
    $hash->{$dinuc}++;
  }

  my $lastnuc = substr( $seq, $hash->{'LENGTH'} - 1, 1 );
  $hash->{$lastnuc}++;

  return $hash;
}
