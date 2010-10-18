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
use File::Temp qw(tempfile);
use Getopt::Long;

###############################################################################
# PROGRAM PATHS
###############################################################################

my $RNAfold        = 'RNAfold';
my $ushuffle       = 'ushuffle';
my $RNAfoldoptions = '-d2 -noPS';

#my $RNAfoldoptions = '-d2 --noPS'; # use for Vienna RNA package 2.0

###############################################################################
# OPTIONS PROCESSING
###############################################################################

my $n = 1000;
my $r = 1;
my $mono;
my $di;
my $help;
my $k;
my $verbose;
my $nochecking;
my $flag       = '';
my @global_avg = ();
my @global_sd  = ();
my @seeds      = ();
my $seq;
my $pipe;
my $ensemble;

GetOptions(
  "n:i"            => \$n,
  "r:i"            => \$r,
  "runs:i"         => \$r,
  "mononucleotide" => \$mono,
  "m"              => \$mono,
  "dinucelotide"   => \$di,
  "d"              => \$di,
  "V"              => \$verbose,
  "verbose"        => \$verbose,
  "help"           => \$help,
  "h"              => \$help,
  "nochecking"     => \$nochecking,
  "pipe"           => \$pipe,
  "s:s"            => \$seq,
  "sequence:s"     => \$seq,
  "ensemble"       => \$ensemble
);

if ( !$nochecking ) {
  my $RNAfold_check  = `which $RNAfold`;
  my $RNAfold_flag   = ( $RNAfold_check or -e $RNAfold ) ? 1 : 0;
  my $ushuffle_check = `which $ushuffle`;
  my $ushuffle_flag  = ( $ushuffle_check or -e $ushuffle ) ? 1 : 0;

  if ( $RNAfold_flag == 0 or $ushuffle_flag == 0 ) {
    print STDERR "ERROR: RNAfold executable not found.\n"  if ( $RNAfold_flag == 0 );
    print STDERR "ERROR: ushuffle executable not found.\n" if ( $ushuffle_flag == 0 );
    exit 1;
  }
}

if ( $mono and $di ) {
  print STDERR "ERROR: These options cannot be combined.\n";
  exit 1;
}

if ( $mono and !$di ) {
  $k = 1;
} elsif ( !$mono and $di ) {
  $k = 2;
} elsif ( !$mono and !$di ) {
  $k = 1;
}

if ( $k < 1 or $k > 2 ) {
  print STDERR "ERROR: only mononucleotide and dinucleotide shuffling allowed.\n";
  exit 1;
}

$flag = ( $k == 2 ) ? 'dinucelotide' : 'mononucleotide';
my $flag_energy = ($ensemble) ? 'ensemble' : 'mfe';
if ($pipe) {
  $seq = <STDIN>;
}

if ( !$seq ) {
  print "ERROR: No sequence provided.\n";
  print usage();
  exit 1;
}

if ($ensemble) {
  $RNAfoldoptions .= ' -p';
}

###############################################################################
# MAIN
###############################################################################

chomp $seq;
if ( $seq =~ m/[^A-Z]/i ) {
  print STDERR "ERROR: Strange sequence $seq.\n";
  exit 1;
}

my $mfe;
my $ensemble_energy;

# call RNAfold to get MFE or ensemble free energy
my @RNAfoldout = `echo $seq | $RNAfold $RNAfoldoptions`;
foreach my $line (@RNAfoldout) {
  if ( $line =~ m/.*\s\(\s*(-?\d+\.\d+)\)\n/ ) {
    $mfe = $1;
  }
  if ( $line =~ m/.*\s\[\s*(-?\d+\.\d+)\]\n/ ) {
    $ensemble_energy = $1;
  }
}

# shuffle 'n fold
for ( 1 .. $r ) {

  # choose a seed between 1 ... 100,000,000
  my $seed = int( rand(100000000) );
  my @out  = `$ushuffle -s $seq -k $k -n $n -seed $seed | $RNAfold $RNAfoldoptions`;

  my $avg    = 0;
  my @values = ();
  foreach my $i ( 0 .. $#out ) {
    if ( $out[$i] =~ m/.*\s\(\s*(-?\d+\.\d+)\)\n/ and !$ensemble ) {
      push @values, $1;
      $avg += $1;
    }
    if ( $out[$i] =~ m/.*\s\[\s*(-?\d+\.\d+)\]\n/ and $ensemble ) {
      push @values, $1;
      $avg += $1;
    }
  }

  if ( $#values != $n - 1 ) {
    print STDERR "ERROR: Failed to generate $n shuffled sequences\n";
    exit 1;
  }

  $avg = $avg / $n;
  my $stdv = 0;
  foreach my $i ( 0 .. $#values ) {
    $stdv += ( $values[$i] - $avg )**2;
  }
  $stdv = $stdv / ( $n - 1 );
  $stdv = sqrt($stdv);
  $stdv = sprintf( "%.5f", $stdv );
  $avg  = sprintf( "%.5f", $avg );
  push @global_avg, $avg;
  push @global_sd,  $stdv;
  push @seeds,      $seed;

}

# calculate z-score
my $z = 0;
my $energy = ($ensemble) ? $ensemble_energy : $mfe;
if ( $r == 1 ) {
  $z = ( $energy - $global_avg[0] ) / $global_sd[0];
} else {
  for my $i ( 0 .. $r - 1 ) {
    $z += ( $energy - $global_avg[$i] ) / $global_sd[$i];
  }
  $z = $z / $r;
}

# output
$z = sprintf( "%.2f", $z );
if ($verbose) {
  print "$z,$flag,$flag_energy,$n,$energy,", join( ";", @seeds ), ",", join( ";", @global_avg ),
    ",",
    join( ";", @global_sd ), ",$seq\n";
} else {
  print "$z\n";
}

sub usage {
  print <<EOF;

NAME

  RNAfoldzscore.pl - Calculates the z-score of a given DNA/RNA sequence using
  RNAfold and ushuffle.


DESCRIPTION

  RNAfoldzscore.pl takes a DNA/RNA sequence as input and calculates a
  z-score as a measure of thermodynamic stability.

  Frist RNAfold is called to caluclate the minimum free energy of the
  native sequence. 'n' shuffled sequences are then generated using
  ushuffle using a mononucleotide or dinucleotide preserving algorithm
  and subsequently folded with RNAfold.

USAGE

   echo ACGTACGTACGT | RNAfoldzscore.pl --pipe [OPTIONS]

   RNAfoldzscore.pl -s ACGTACGTACGT [OPTIONS]


OPTIONS

  -n                     Number of shuffled sequences to be generated.
                         Default: 1000

  -r, --runs             Number of runs to use. Z-scores will be averaged
                         r times. Default: 1

  -m, --mononucleotide   Mononucelotide background model.
  
  -d, --dinucleotide     Dinucelotide background model.

  --nochecking           Do not check if RNAfold and ushuffle executables
                         are found.

  --pipe                 RNA/DNA sequence is read from STDIN instead of
                         the paramter '-s' or '--sequence'.

  -s, --sequence         RNA/DNA sequence.

  -V, --verbose          Gives verbose output of the form:
                         z-score,model,n,mfe,seed(s),averages,stdvs

  -h, --help             Display this help message.

  --ensemble             Calculate the z-score based on the ensemble free
                         energy instead of the minimum free energy.

EOF
  exit;
}
