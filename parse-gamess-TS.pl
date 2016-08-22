#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
#use Excel::Writer::XLSX;
#use Spreadsheet::WriteExcel;
#use Chemistry::OpenBabel;

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


my @energy;
my @systemname;
my @imaginarymode;

FN1:
while ($fn = <$filelist>) {
    chomp $fn;

    my $filedat="$fn.dat";
    open my $Rd, '<', $filedat or goto FN1;
    print "parsing $filedat\n";
    push @systemname, $fn;

    my $tmpen;
    my $tmpmode;
    my $tmp1; my $tmp2;
    while (<$Rd>) {
	chomp;
	if (/^ENERGY IS/) {
	    my @tmparray = split(/\s+/);
	    $tmpen = $tmparray[2];
	}
	if (/^MODE    1 /) {
	    $_ =~ s/MODE    1   FREQUENCY=//;
	    ($tmpmode,$tmp1) = split(/\s+/);
	}
    }

    push @energy, $tmpen;
    push @imaginarymode, $tmpmode;

    close $Rd;
}


seek $filelist, 0, 0;

my @maxgrad;
my @rmsgrad;

FN2:
while ($fn = <$filelist>) {
    chomp $fn;
    my $fileout="$fn.out";
    open my $R, '<', $fileout or goto FN2;
    print "parsing $fileout\n";

    my @tmpgrad;
    while(<$R>) {
	if (/MAXIMUM GRADIENT =/){
	    chomp;
	    @tmpgrad = split(/\s+/);
	}
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

    push @maxgrad,$tmpgrad[4];
    push @rmsgrad,$tmpgrad[8];

    close $R;

    my $writefile="$fn.xyz";
    open my $W, '>', $writefile or die "$writefile: $!\n";

    print $W "$i\n";
    print $W "$fileout: $searchline";
    for(my $k=0; $k < $i; $k++) {
	my ($tmp1,$tmp2,$tmp3,$tmp4,$tmp5) = split(/\s+/, $line[$k]);
        print $W "$tmp1  $tmp3  $tmp4  $tmp5\n";
    }
    print $W "\n";
    close $W;

}


#print @systemname, "\n";
#print @energy, "\n";
#print @imaginarymode,"\n";
#print @maxgrad,"\n";
#print @rmsgrad,"\n";

print "\n\n****\n";
my $firstline = "$list_file_name";
$firstline =~ s/_/ /g;
print "$firstline\n";
my $hostname = "";
   $hostname .= hostname;
   $hostname .= ":";
   $hostname .= getcwd;
print $hostname,"\n\n";

print "System name\t";
print "Transition Energy\t";
print "Imaginary mode (cm^-1)\t";
print "Maximum gradient\t";
print "RMS gradient\n";
my $jj=0;
foreach (@systemname) {
    print $systemname[$jj], "\t";
    print $energy[$jj], "\t";
    print $imaginarymode[$jj], "\t";
    print $maxgrad[$jj], "\t";
    print $rmsgrad[$jj], "\t";
    print "\n";
    $jj++;
}
print "\n";

open my $W, '>', "$list_file_name.csv" or die "$list_file_name.csv: $!\n";

print $W "$firstline\n";

print $W $hostname,"\n";
print $W "System name,";
print $W "Transition Energy,";
print $W "Imaginary mode (cm^-1),";
print $W "Maximum gradient,";
print $W "RMS gradient\n";
$jj=0;
foreach (@systemname) {
    print $W $systemname[$jj], ",";
    print $W $energy[$jj], ",";
    print $W $imaginarymode[$jj], ",";
    print $W $maxgrad[$jj], ",";
    print $W $rmsgrad[$jj];
    print $W "\n";
    $jj++;
}

close $W;

#my $workbook = Excel:Writer::XLSX->new("$list_file_name.xlsx");
#my $worksheet = $workbook->add_worksheet();

#$worksheet->write(1,0,

