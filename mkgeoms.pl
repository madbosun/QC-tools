#!/usr/bin/perl
use strict;
use warnings;

if (!$ARGV[0]) {
    print "./mkgeoms.pl filelist ext#\n";
    exit 1;
}

my $list_file_name=$ARGV[0];
my $extnum=$ARGV[1];
my $fn="";

my @sym = qw(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti Y Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn);
my $searchline="";
my $i=0;
my @line;
my $group="ccc";

open my $filelist, '<', $list_file_name or die "$list_file_name: $!\n";

FN:
while ($fn = <$filelist>) {
    chomp $fn;
    my $file="$fn.$extnum.out";
    open my $R, '<', $file or goto FN;
    print "parsing $file\n";

    while(<$R>) {
	if (/THE POINT GROUP OF THE/) {
	    chomp;
	    $_ =~ s/ THE POINT GROUP OF THE MOLECULE IS //;
	    $group=$_;
	}
	if (/^ NSERCH/) {
	    $_ =~ s/=\s+/=/g;
	    $searchline = $_;
	    }
	if (/COORDINATES OF ALL ATOMS ARE/) {
	    $_ = <$R>; $_ = <$R>;
	    $i = 0;
	    while(<$R>) {
		last if(/^\n/);
		$_ =~ s/^ //;
		$_ =~ s/        //;
		$line[$i++] = $_;
	    }
	}
    }

    close $R;

    my $writefile="$fn.gamin";
    open my $W, '>', $writefile or die "$writefile: $!\n";

    print $W " \$data\n";
    print $W "$file: $searchline";
    print $W " $group\n";
    for(my $k=0; $k < $i; $k++) {
    print $W $line[$k]
    }
    print $W " \$end\n\n";
    close $W;

}
