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
#
# computes vibrational mode based on calculations coming out of gms-displace.pl
#
#################################################################

print STDERR "\n-- BEWARE: this script assumes that all outputs were created from gms-displace.pl\n";
print STDERR   "           You are taking chances in life if you hand created in and output files\n\n";

my $deltas= 0.0;

use Getopt::Long;
my $filename=$ARGV[0];

if (@$ARGV[0]) {
  print "WARNING: no input file name.\n";
  exit 1;
}

my $moleculename=$ARGV[1];

if (!$moleculename) {
  print "WARNING: no base molecule name.\n";
  exit 1;
}

(my $rootname = $filename) =~ s/\.[^.]+$//; $rootname =~ s{.*/}{};

if ($filename =~ /.gz$/) {
    $rootname =~ s/\.[^.]+$//;
    open R, '-|', 'gunzip', '-c', $filename 
        or die "Could not read from $filename: $!";
} else {
    open(R,"<$filename") or die "cannot open $filename\n";
}

my $reffreq;
my $redmass;
my $amumass;
my @inertia;
my @mp2mode;

my @sym = qw(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar);
my @atmaslist = qw(1.0078250 4.002602 6.941 9.012182 10.811 12.0107 14.0067 15.9994 18.9984032 20.1797 22.989770 24.3050 26.981538 28.0855 30.973761 32.065 35.453 39.948);
my %masshash;
for(my $k=0; $k<scalar(@sym); $k++) {
    $masshash{$sym[$k]} = [] unless exists $masshash{$sym[$k]};
    $masshash{$sym[$k]} = $atmaslist[$k];
}

print "-" x 64,"\n";
print STDERR "  Parsing $filename\n";
print "  Root GAMESS file is $rootname for the molecule $moleculename\n";

FILEREAD:
while(<R>) {
    if (/FREQUENCY/) {
        $_ =~ s/^\s+//;
        my @ftmp = split(/\s+/);
        $reffreq = $ftmp[1];
        my $itmp = $ftmp[2];
        print "  The Reference mode is: $reffreq\n";
# symmetry line
        $_ = <R>;
	    ($_ = <R>) =~ s/^\s+//;
        my @mtmp = split(/\s+/);
	    #$redmass = Math::BigFloat->new($mtmp[2]);
	    $redmass = $mtmp[2];
	    $amumass = $redmass*1822.89;
        print "  The Reduced mass is (AMU): $amumass\n";

        while(<R>) {
	    if (/MODE FREQ/) {
		my $i = 0;
		while(<R>) {
		    last if (/^\n/);
		    chomp;
		    $_ =~ s/^\s+//;
		    $mp2mode[$i++]=$_;
		}
	    }
	    if (/THE MOMENTS/){
		($_ = <R>) =~ s/^\s+//;
		@inertia = split(/\s+/);
		printf "  The moments of inertial (AMU*BOHR**2): %10.5f %10.5f %10.5f\n", @inertia;
                last FILEREAD;
	    }
	}

    }
}

open(Rmds,"<$rootname-mds.xyz") or die "cannot open $rootname-mds.xyz\n";

$_ = <Rmds>; ($_ = <Rmds>) =~ s/^\s+//;
my @tmpline = split(/\s+/);
$deltas = $tmpline[4];
print "  The read DS (ang): $deltas\n";
my $delta=$deltas/0.529177249;

close(Rmds);

my $Emds = 0.0;
my $Eref = 0.0;
my $Epds = 0.0;

open(Omds,"<$rootname-mds.out") or die "cannot open $rootname-mds.out\n";
while (<Omds>){ if (/    0    /) { my @tmpline = split(/\s+/); $Emds = $tmpline[3]; } }
close(Omds);
open(Oref,"<$rootname-ref.out") or die "cannot open $rootname-ref.out\n";
while (<Oref>){ if (/    0    /) { my @tmpline = split(/\s+/); $Eref = $tmpline[3]; } }
close(Oref);
open(Opds,"<$rootname-pds.out") or die "cannot open $rootname-pds.out\n";
while (<Opds>){ if (/    0    /) { my @tmpline = split(/\s+/); $Epds = $tmpline[3]; } }
close(Opds);

printf "  The TPF MP2 energies (a.u.): %24.14f %24.14f %24.14f\n", $Emds,$Eref,$Epds;

my $D2E_Ds2 = ($Emds + $Epds - 2.0*$Eref)/($delta*$delta);
my $sign = 1.0;
if ($D2E_Ds2 < 0) { $sign = -1.0; }
my $mode = $sign*sqrt(abs($D2E_Ds2)/$amumass)*219474.63067;
printf "-- CHECK: MP2 TPF and error is (cm**-1): %10.5f, %10.5f\n",$mode,$mode+$reffreq;

$Emds = 0.0;
$Eref = 0.0;
$Epds = 0.0;

open(Omds,"<$rootname-mds.out") or die "cannot open $rootname-mds.out\n";
while (<Omds>){ if (/Final CCSD/) { my @tmpline = split(/\s+/); $Emds = $tmpline[5]; } }
close(Omds);
open(Oref,"<$rootname-ref.out") or die "cannot open $rootname-ref.out\n";
while (<Oref>){ if (/Final CCSD/) { my @tmpline = split(/\s+/); $Eref = $tmpline[5]; } }
close(Oref);
open(Opds,"<$rootname-pds.out") or die "cannot open $rootname-pds.out\n";
while (<Opds>){ if (/Final CCSD/) { my @tmpline = split(/\s+/); $Epds = $tmpline[5]; } }
close(Opds);

printf "\n  The TPF CCSD energies (a.u.): %24.14f %24.14f %24.14f\n", $Emds,$Eref,$Epds;

$D2E_Ds2 = ($Emds + $Epds - 2.0*$Eref)/($delta*$delta);
$sign = 1.0;
if ($D2E_Ds2 < 0) { $sign = -1.0; }
$mode = $sign*sqrt(abs($D2E_Ds2)/$amumass)*219474.63067;
printf "  The TPF (cm**-1): %10.5f\n\n",$mode;

$rootname =~ s/g2cc_tpf_//;

open(WSAM,">$moleculename-$rootname.sam3") or die "cannot open $moleculename-$rootname.sam3: $!\n";

print WSAM "$moleculename $rootname\n";
print WSAM "-" x 64,"\n\n";
print WSAM  "Multiplicity = 0\n";
print WSAM "Barrier = __ (kcal/mol)\n";
printf WSAM "MOI XX=%11.5f YY=%11.5f ZZ=%11.5f (AMU-Bohr**2)\n", @inertia;
printf WSAM "TPF = %10.5f (cm**-1)\n\n", $mode;
print WSAM "Harmonic modes for $moleculename, TS $rootname using MP2/6-31++G** @ LCCD/6-31++G**\n","-" x 64,"\n";

my @tmp = split /\s+/, $mp2mode[0];
printf WSAM " %3u %10.3f\n", $tmp[0], -1*$tmp[1],"\n";
for (my $i = 1; $i < scalar(@mp2mode); ++$i) {
    my @tmp = split /\s+/, $mp2mode[$i];
    printf WSAM " %3u %10.3f\n", $tmp[0], $tmp[1],"\n";
}

close WSAM;

exit 0;


