#!/usr/bin/perl -w
use strict;
#use Math::BigFloat;

#################################################################
#
# Author: Jason Byrd
# Contact: byrd.jason@ensco.com
# date: Fall 2016
#
#################################################################
# read a gamess hessian output file, generate cartesian coordinates of
# displacements along the normal mode vector
#################################################################

my $deltas= 3.0; #Math::BigFloat->new(3.0);

use Getopt::Long;
GetOptions("ds:s" => \$deltas);
my $filename=$ARGV[0];

if (@$ARGV[0]) {
  print "WARNING: no input file name.\n";
  exit 1;
}

$deltas *= 0.0001;

print "  Using a delta s of $deltas\n";

(my $rootname = $filename) =~ s/\.[^.]+$//; $rootname =~ s{.*/}{};

open(R,"<$filename") or die "cannot open $filename\n";

my @at;
my @x;
my @y;
my @z;

my @XX;
my @YY;
my @ZZ;

my $reffreq;
my $redmass;

my @sym = qw(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar);
my @atmaslist = qw(1.0078250 4.002602 6.941 9.012182 10.811 12.0107 14.0067 15.9994 18.9984032 20.1797 22.989770 24.3050 26.981538 28.0855 30.973761 32.065 35.453 39.948);
my %masshash;
for(my $k=0; $k<scalar(@sym); $k++) {
    $masshash{$sym[$k]} = [] unless exists $masshash{$sym[$k]};
    $masshash{$sym[$k]} = $atmaslist[$k];
}

FILEREAD:
while(<R>) {
# get standard coordinates from gamess output
	#if (/COORDINATES OF ALL ATOMS ARE/) {
	if (/ ATOM      ATOMIC                      COORDINATES/) {
	    $_ = <R>; #$_ = <R>;
	    my $i = 0;
	    while(<R>) {
            last if(/^\n/);
            last if(/^\s+[^\w]/);
            $_ =~ s/^ //;
            $_ =~ s/        //;
            my @tmp  = split(/\s+/);
            $at[$i] = $tmp[0];
            #$x[$i] = Math::BigFloat->new($tmp[2]);
            #$y[$i] = Math::BigFloat->new($tmp[3]);
            #$z[$i] = Math::BigFloat->new($tmp[4]);
            $x[$i] = $tmp[2]*0.529177249;
            $y[$i] = $tmp[3]*0.529177249;
            $z[$i] = $tmp[4]*0.529177249;
            ++$i;
	    }
	}
    if (/FREQUENCY/) {
        $_ =~ s/^\s+//;
        my @ftmp = split(/\s+/);
        $reffreq = $ftmp[1];
        my $itmp = $ftmp[2];
        print "  The Reference mode is: $reffreq\n";
        if ($itmp =~ /I/) {
            print "    -- Coordinates are confirmed to be at least an order 1 saddle\n";
        } else {
            print "    -- Coordinates are NOT at a saddle\n";
            exit 1;
        }
# symmetry line
        $_ = <R>;
	    ($_ = <R>) =~ s/^\s+//;
        my @mtmp = split(/\s+/);
	    #$redmass = Math::BigFloat->new($mtmp[2]);
	    $redmass = $mtmp[2];
        print "  The Reduced mass is: $redmass\n";

# ir intensity
        $_ = <R>;
# blank line
        $_ = <R>;
        my $i = 0;

        while (<R>) {
            #chomp;
            if(/^\s*\n/) {last FILEREAD;}
            $_ =~ s/^\s+//;

            my @tmp1 = split(/\s+/);
            #$XX[$i] = Math::BigFloat->new($tmp1[3]);
            $XX[$i] = $tmp1[3];

            ($_ = <R>) =~ s/^\s+//;
            my @tmp2 = split(/\s+/);
            #$YY[$i] = Math::BigFloat->new($tmp2[2]);
            $YY[$i] = $tmp2[2];

            ($_ = <R>) =~ s/^\s+//;
            my @tmp3 = split(/\s+/);
            #$ZZ[$i] = Math::BigFloat->new($tmp3[2]);
            $ZZ[$i] = $tmp3[2];

       #     print "$XX[$i] $YY[$i] $ZZ[$i]\n";
            ++$i;
        }

# this is where we finish assuming first mode only
        last FILEREAD;
    }
}


print "  Assuming standard isotope mass\n";
my $totalmass=0.0;
my $cmX=0.0;
my $cmY=0.0;
my $cmZ=0.0;
for (my $k=0; $k < scalar(@at); ++$k) {
    $totalmass += $masshash{$at[$k]};
    $cmX += $masshash{$at[$k]}*$x[$k];
    $cmY += $masshash{$at[$k]}*$y[$k];
    $cmZ += $masshash{$at[$k]}*$z[$k];
}

$cmX /= $totalmass;
$cmY /= $totalmass;
$cmZ /= $totalmass;

print "\n  Normal Mode 1 displacement vector\n\n";
print "   At       H_X_           H_Y_           H_Z_\n";
for (my $k=0; $k < scalar(@at); ++$k) {
    printf "    %-4s %15.10f %15.10f %15.10f\n", $at[$k],$XX[$k],$YY[$k],$ZZ[$k];
}
print "\n";

#my $norm = Math::BigFloat->bzero();
my $norm = 0.0;
for (my $k=0; $k < scalar(@at); ++$k) {
    $norm += $XX[$k]*$XX[$k] + $YY[$k]*$YY[$k] + $ZZ[$k]*$ZZ[$k];
}
my $delta = sqrt($deltas*$deltas/$norm);

my @Px;
my @Py;
my @Pz;
my @Mx;
my @My;
my @Mz;
for (my $k=0; $k < scalar(@at); ++$k) {
# compute +ds
    $Px[$k] = $x[$k] + $XX[$k]*$delta;
    $Py[$k] = $y[$k] + $YY[$k]*$delta;
    $Pz[$k] = $z[$k] + $ZZ[$k]*$delta;

# compute -ds
    $Mx[$k] = $x[$k] - $XX[$k]*$delta;
    $My[$k] = $y[$k] - $YY[$k]*$delta;
    $Mz[$k] = $z[$k] - $ZZ[$k]*$delta;

# shift to TS center of mass    
    $x[$k] -= $cmX;
    $y[$k] -= $cmY;
    $z[$k] -= $cmZ;
    $Px[$k] -= $cmX;
    $Py[$k] -= $cmY;
    $Pz[$k] -= $cmZ;
    $Mx[$k] -= $cmX;
    $My[$k] -= $cmY;
    $Mz[$k] -= $cmZ;
}

open(Wp,">$rootname-pds.xyz");
open(Wr,">$rootname-ref.xyz");
open(Wm,">$rootname-mds.xyz");

print Wp scalar(@at),"\n";
print Wr scalar(@at),"\n";
print Wm scalar(@at),"\n";

print Wp "$rootname + ds = $deltas reduced mass = $redmass\n";
print Wr "$rootname reference reduced mass = $redmass\n";
print Wm "$rootname - ds = $deltas reduced mass = $redmass\n";

print "\n  Reference cartesian coordinates\n\n";
print "   At       _X_             _Y_             _Z_\n";
for (my $k=0; $k < scalar(@at); ++$k) {
    printf "    %-4s %15.10f %15.10f %15.10f\n", $at[$k],$x[$k],$y[$k],$z[$k];
    printf Wr "%-4s %15.10f %15.10f %15.10f\n", $at[$k],$x[$k],$y[$k],$z[$k];
}
print "\n";

print "\n  Plus ds cartesian coordinates\n\n";
print "   At       _X_             _Y_             _Z_\n";
for (my $k=0; $k < scalar(@at); ++$k) {
    printf "    %-4s %15.10f %15.10f %15.10f\n", $at[$k],$Px[$k],$Py[$k],$Pz[$k];
    printf Wp "%-4s %15.10f %15.10f %15.10f\n", $at[$k],$Px[$k],$Py[$k],$Pz[$k];
}
print "\n";

print "\n  Minus ds cartesian coordinates\n\n";
print "   At       _X_             _Y_             _Z_\n";
for (my $k=0; $k < scalar(@at); ++$k) {
    printf "    %-4s %15.10f %15.10f %15.10f\n", $at[$k],$Mx[$k],$My[$k],$Mz[$k];
    printf Wm "%-4s %15.10f %15.10f %15.10f\n", $at[$k],$Mx[$k],$My[$k],$Mz[$k];
}
print "\n";

print Wp "\n";
print Wr "\n";
print Wm "\n";

close(Wp);
close(Wr);
close(Wm);

my $deltasau = $deltas/(0.52917706*0.52917706);

