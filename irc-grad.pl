#!/usr/bin/perl -w
use strict;

if (@$ARGV[0]) {
  print "WARNING: no input file name.\n";
  exit 1;
}

my @sym = qw(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti Cr Mn Fe Co Ni V Cu Zn Ga Ge As Se Br Kr Rb);
my $k;
my %symhash;
for($k=1; $k<scalar(@sym); $k++) {
          $symhash{$sym[$k]} = [] unless exists $symhash{$sym[$k]};
          $symhash{$sym[$k]} = $k;
}

my $filename=$ARGV[0];

open(Ro,"$filename.out") || die "cannot open $filename.out\n";
open(Rd,"$filename.dat") || die "cannot open $filename.dat\n";


# getting DATA and GRAD from .out file, lst one per POINT is the required block
while(<Ro>) {
    my $prefix=0;
    if (/CONSTRAINED OPTIMIZATION POINT/) {
        my@temp = split(/\s+/);
        $prefix = $temp[8];
        print "$filename.$prefix\n";
        open(W,">$filename.$prefix");
        print "$filename.$prefix\n";

        print W '$DATA',"\n";
        print W "$filename IRC POINT $prefix\n";
        $_ = <Ro>; $_ = <Ro>; $_ = <Ro>;
        while(<Ro>) {
            if(/^ \.\./){
                print W '$END',"\n";
                last;
            }
            $_ =~ s/^ //;
            print W $_;
        }
    }
    if (/ UNITS ARE HARTREE/) {
        print W '$GRAD',"\n";
        print W "$filename IRC POINT $prefix\n";

      my $x;
      my $y;
      my $z;
      my $at;
      $k=1;

        while(<Ro>) {
            if(/^\n/) {
                print W '$END',"\n";
                close(W);
                last;
            }
            $_ =~ s/^\s+\d+\s+//;
            ($at,$x,$y,$z) = split;
            printf W "%-4s %2s. %12s %12s %12s\n", $at, $symhash{$at}, $x, $y, $z;
            $k++;
        }
    }
}
my $count=0;
while(<Rd>) {
    if (/\$HESS/) {
        $count++;
        open(W,">>$filename.$count");
        while (<Rd>){
            print W $_;
            if (/\$END/){
                close W;
                last;
            }
        }
    }
}

