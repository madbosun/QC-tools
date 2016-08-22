#!/usr/bin/perl
#

  use strict;
  use Getopt::Long qw(:config bundling);
  use warnings qw(void syntax severe redefine);
  use warnings FATAL => qw(io overflow);
  no warnings qw(uninitialized);

  my $xyz="";
  my $zxy="";
  my $yzx="";
  my $help="";

  my @sym = qw(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti Cr Mn Fe Co Ni V Cu Zn Ga Ge As Se Br Kr Rb);
  my $k;
  my %symhash;
  for($k=1; $k<scalar(@sym); $k++) {
          $symhash{$sym[$k]} = [] unless exists $symhash{$sym[$k]};
          $symhash{$sym[$k]} = $k;
  }

  my $x;
  my $y;
  my $z;
  my $at;
  $k=1;
  while(<>) {
          chomp;
        ($at,$x,$y,$z) = split;
            printf "%-4s %2s %12s %12s %12s\n", $at, $symhash{$at}, $x, $y, $z;
          $k++;
  }

