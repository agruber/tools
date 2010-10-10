#!/usr/bin/perl

# Copyright 2009, 2010 Andreas R. Gruber
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.  You should have received a
# copy of the GNU General Public License along with Paperpile.  If
# not, see http://www.gnu.org/licenses.

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile);
use Term::ANSIColor;

###############################################################################
# PROGRAM PATHS AND SANITY CHECKS
###############################################################################

my $clustalw       = 'clustalw';
my $clustalw_check = `which $clustalw`;
my $clustalw_flag  = ( $clustalw_check or -e $clustalw ) ? 1 : 0;
my $dialign        = 'dialign2-2';
my $dialign_check  = `which $dialign`;
my $dialign_flag   = ( $dialign_check or -e $dialign ) ? 1 : 0;
my $mafft          = 'mafft';
my $mafft_check    = `which $mafft`;
my $mafft_flag     = ( $mafft_check or -e $mafft ) ? 1 : 0;
my $locarna        = 'mlocarna';
my $locarna_check  = `which $locarna`;
my $locarna_flag   = ( $locarna_check or -e $locarna ) ? 1 : 0;

print STDERR "========== Program stettings ==========\n";
if ( $clustalw_flag == 1 ) {
  print STDERR "Clustalw (";
  print STDERR color 'bold green';
  print STDERR "active";
  print STDERR color 'reset';
  print STDERR "): $clustalw";
  print STDERR color 'bold green';
  print STDERR "\tLETTER: C\n";
  print STDERR color 'reset';
} else {
  print STDERR "Clustalw (";
  print STDERR color 'bold red';
  print STDERR "inactive";
  print STDERR color 'reset';
  print STDERR "): $clustalw\n";
}
if ( $dialign_flag == 1 ) {
  print STDERR "DIALIGN  (";
  print STDERR color 'bold green';
  print STDERR "active";
  print STDERR color 'reset';
  print STDERR "): $dialign";
  print STDERR color 'bold green';
  print STDERR "\tLETTER: D\n";
  print STDERR color 'reset';
} else {
  print STDERR "DIALIGN  (";
  print STDERR color 'bold red';
  print STDERR "inactive";
  print STDERR color 'reset';
  print STDERR "): $dialign\n";
}
if ( $mafft_flag == 1 ) {
  print STDERR "MAFFT    (";
  print STDERR color 'bold green';
  print STDERR "active";
  print STDERR color 'reset';
  print STDERR "): $mafft";
  print STDERR color 'bold green';
  print STDERR "\tLETTER: M\n";
  print STDERR color 'reset';
} else {
  print STDERR "MAFFT    (";
  print STDERR color 'bold red';
  print STDERR "inactive";
  print STDERR color 'reset';
  print STDERR "): $mafft\n";
}
if ( $locarna_flag == 1 ) {
  print STDERR "LocARA   (";
  print STDERR color 'bold green';
  print STDERR "active";
  print STDERR color 'reset';
  print STDERR "): $locarna";
  print STDERR color 'bold green';
  print STDERR "\tLETTER: L\n";
  print STDERR color 'reset';
} else {
  print STDERR "LocARNA  (";
  print STDERR color 'bold red';
  print STDERR "inactive";
  print STDERR color 'reset';
  print STDERR "): $locarna\n";
}

###############################################################################
# MAIN
###############################################################################

my $input_stk     = $ARGV[$#ARGV];
my $current_file  = $input_stk;
my $iteration     = 0;
my @USED_FILENAES = ();
my $letter        = 'X';
my @orig_seqs     = ();
my $help;

GetOptions(
  "letter:s" => \$letter,
  "l:s"      => \$letter,
  "help"     => \$help,
  "h"        => \$help
);

if ( $letter !~ m/^[A-Z]$/ ) {
  print STDERR "Letter symbol '$letter' is not supported to indicate ";
  print STDERR "the stretch of nucleotides to be realigned.\n";
  print STDERR "Do nothing and exit.\n";
  exit 1;
}

if ( ! $input_stk or $help ) {
  usage();
  exit 1;
}


if ( !-e $input_stk ) {
  usage();
  exit 1;
}

print STDERR "Using symbol '$letter' to find stretches of nucleotides ";
print STDERR "that should be realigned.\n\n";

while (1) {
  $iteration++;

  my @seqs  = ();
  my $start = -1;
  my $end   = -1;

  ############################################################################
  # read in the stockholm file and scan the "GC Anno" line for the next
  # stretch of $letter. Start and end of the $letter stretch will be stored.
  ############################################################################

  open( INFILE, $current_file );
  while ( my $line = <INFILE> ) {
    chomp $line;
    next if ( $line =~ m/STOCKHOLM/ );
    next if ( $line !~ m/[A-Z]/i );
    if ( $line =~ m/^(#=GC\sAnno)\s+(\S+)/ ) {
      my $tmp = {};
      $tmp->{name} = $1;
      $tmp->{seq}  = $2;
      push @seqs, $tmp;
      my @temp = split( //, $tmp->{seq} );

      # search for $letter stretch
      for ( my $i = 0 ; $i <= $#temp ; $i++ ) {
        $start = $i if ( $temp[$i] eq $letter and $start == -1 );
        $end = $i if ( $temp[$i] ne $letter and $start > -1 and $end == -1 );
      }
      if ( $start > -1 and $end == -1 and $temp[$#temp] eq $letter ) {
        $end = $#temp;
      }

    } elsif ( $line !~ m/^#/ ) {
      my $tmp = {};
      ( my $name, my $seq ) = split( /\s+/, $line );
      $tmp->{name} = $name;
      $tmp->{seq}  = $seq;
      push @seqs, $tmp;
    } elsif ( $line =~ m/^(#=GC.*)\s+(\S+)/ ) {
      my $tmp = {};
      $tmp->{name} = $1;
      $tmp->{seq}  = $2;
      push @seqs, $tmp;
    }
  }
  close(INFILE);

  # check if sequences are still the same
  if ( $iteration == 1 ) {
    @orig_seqs = @seqs;
  } else {
    for my $i ( 0 .. $#orig_seqs ) {
      next if ( $orig_seqs[$i]->{name} =~ m/^#/ );
      ( my $seqA = $orig_seqs[$i]->{seq} ) =~ s/(-|\.)//g;
      ( my $seqB = $seqs[$i]->{seq} )      =~ s/(-|\.)//g;
      if ( uc($seqA) ne uc($seqB) ) {
        print STDERR "There seems to be an inconsistancy.\n";
        print STDERR "Sequence: ", $orig_seqs[$i]->{name}, "\n";
        print STDERR "$seqA\n$seqB\n";
        print STDERR "Do nothing and exit.\n\n";
        exit 1;
      }
    }
  }

  ############################################################################
  # If we cannot find a X stretch $start will be -1. In this case we print the
  # alignment that is stored in $current_file and exit.
  ############################################################################

  if ( $start == -1 ) {

    open( FILE_FINAL, $current_file );
    while ( my $line = <FILE_FINAL> ) {
      print $line;
    }

    # cleaning temp files
    foreach my $file (@USED_FILENAES) {
      `rm $file`;
      `rm $file.aln`    if ( -e "$file.aln" );
      `rm -r $file.out` if ( -e "$file.out/results/result.tree" );
    }

    exit 0;
  }
  print STDERR "Found a stretch to realign from columns $start - $end.\n";

  ############################################################################
  # Create a temporary  FASTA styled file as input for the alignment
  # programs.
  ############################################################################

  ( my $fh, my $tmpfile ) = tempfile();
  push @USED_FILENAES, $tmpfile;
  my $max = 0;
  my $min = 10e6;
  foreach my $seq (@seqs) {
    $max = length( $seq->{name} ) if ( length( $seq->{name} ) > $max );
    next if ( $seq->{name} =~ m/^#/ );
    my @tmp = split( //, $seq->{seq} );
    my $part = '';
    for ( my $i = $start ; $i < $end ; $i++ ) {
      $part .= $tmp[$i];
    }
    $part =~ s/-//g;
    $part =~ s/\.//g;
    $min = length($part) if ( length($part) < $min );
    if ( length($part) > 0 ) {
      print $fh ">", $seq->{name}, "\n";
      print $fh "$part\n";
    }
  }
  close($fh);

  ############################################################################
  # Call the alignment programs.
  ############################################################################

  my @options = ();
  push @options, 'C' if ( $clustalw_flag == 1 );
  push @options, 'D' if ( $dialign_flag == 1 );
  push @options, 'M' if ( $mafft_flag == 1 and $min > 5 );
  push @options, 'L' if ( $locarna_flag == 1 );

  my $choice = '';
  while ( $choice eq '' ) {
    print STDERR "\tYour choice (", join( "|", @options ), "): ";
    my $answer = <STDIN>;
    chomp $answer;
    foreach my $option (@options) {
      if ( uc($answer) eq $option ) {
        $choice = $option;
        last;
      }
    }
  }

  if ( $choice eq 'C' ) {
    `$clustalw -INFILE=$tmpfile -OUTORDER=INPUT`;
    print STDERR "\tRealigning was done with ClustalW.\n";
  }
  if ( $choice eq 'D' ) {
    `$dialign -cw -n $tmpfile;mv $tmpfile.cw $tmpfile.aln`;
    print STDERR "\tRealigning was done with DIALIGN.\n";
  }
  if ( $choice eq 'M' ) {
    `$mafft --maxiterate 1000 --clustalout $tmpfile > $tmpfile.aln 2>/dev/null`;
    print STDERR "\tRealigning was done with MAFFT.\n";
  }
  if ( $choice eq 'L' ) {
    `$locarna $tmpfile 2>/dev/null`;
    `mv $tmpfile.out/results/result.aln $tmpfile.aln`;
    print STDERR "\tRealigning was done with LocARNA.\n";
  }

  ############################################################################
  # Read in the alignment ouput in CLUSTAL W format.
  ############################################################################

  my %newseqs = ();
  foreach my $entry (@seqs) {
    $newseqs{ $entry->{name} } = '';
  }
  open( FILE_ALN, "$tmpfile.aln" );
  while ( my $line = <FILE_ALN> ) {
    chomp $line;
    next if ( $line !~ m/[A-Z]/i );
    next if ( $line =~ m/DIALIGN/ );
    if ( $line =~ m/(\S+)\s+(\S+)/ ) {
      my $name = $1;
      my $seq  = $2;
      if ( !$newseqs{$name} ) {

        # locarna changes . to _; can cause troubles
        # in recognizng the sequence again
        $name =~ s/_/\./g;
      }
      if ( defined $newseqs{$name} ) {
        $newseqs{$name} .= $seq;
      }
    }
  }

  ############################################################################
  # Get length of the alignment, and some sanity check.
  ############################################################################
  my $length = 0;
  my @l      = ();
  foreach my $key ( keys %newseqs ) {
    $length = length( $newseqs{$key} );
    foreach my $entry (@l) {
      if ( $entry != $length ) {
        print STDERR "Alignment could not be read correctly.\n";
        print STDERR "Procedure failed.\n";
        exit 1;
      }
    }
  }

  ############################################################################
  # Generate the new alignment and write it to a file.
  ############################################################################

  ( my $fh_new, my $tmpfile_new ) = tempfile();
  print $fh_new "# STOCKHOLM 1.0\n\n";
  push @USED_FILENAES, $tmpfile_new;
  foreach my $seq (@seqs) {
    my $newseq = '';

    # regular sequence
    if ( $seq->{name} !~ m/^#/ ) {
      my @tmp = split( //, $seq->{seq} );
      for ( my $i = 0 ; $i < $start ; $i++ ) {
        $newseq .= $tmp[$i];
      }
      if ( $newseqs{ $seq->{name} } ) {

        #print STDERR $newseqs{ $seq->{name} }, "\n";
        $newseq .= uc( $newseqs{ $seq->{name} } );
      } else {
        for ( my $i = 0 ; $i < $length ; $i++ ) {
          $newseq .= '-';
        }
      }
      for ( my $i = $end ; $i <= $#tmp ; $i++ ) {
        $newseq .= $tmp[$i];
      }
    } else {
      my @tmp = split( //, $seq->{seq} );
      for ( my $i = 0 ; $i < $start ; $i++ ) {
        $newseq .= $tmp[$i];
      }
      for ( my $i = 0 ; $i < $length ; $i++ ) {
        $newseq .= '.';
      }
      if ( $end == $#tmp ) {
        $newseq .= '.';
      } else {
        for ( my $i = $end ; $i <= $#tmp ; $i++ ) {
          $newseq .= $tmp[$i];
        }
      }
    }
    print $fh_new $seq->{name};
    for ( my $i = length( $seq->{name} ) ; $i <= $max ; $i++ ) {
      print $fh_new " ";
    }
    print $fh_new $newseq, "\n";
  }
  $current_file = $tmpfile_new;
}

sub usage {
  print <<EOF;

NAME

  RealignStockholm.pl - Takes an alignment in STOCKHOLM format and realigns
  parts that have been marked with 'X' in a comment line called GC Anno.


DESCRIPTION

  RealignStockholm.pl takes an alignment in STOCKHOLM format as outlined
  below and realignes regions marked with X with the alignment program of choice.

  # STOCKHOLM 1.0

  species1     AAACCC---TAGTAG--GGGTTT
  species2     AAACCCTAG--TAG---GGGTTT
  species3     AAACCCT---AG--TAGGGGTTT
  #=GC SS_cons <<<<<<...........>>>>>>
  #=GC Anno    ......XXXXXXXXXXX......

  If alignment programs are not in the path, please edit lines 26 - 37 in
  this Perl script.

USAGE

   RealignStockholm.pl [OPTIONS] input.stk > output.stk

OPTIONS

  -l, --letter           Letter that identifies regions that are
                         to be realigned. Default: X
 
  -h, --help             Display this help message.

  
EOF
  exit;
}
