#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
#use Excel::Writer::XLSX;
#use Spreadsheet::WriteExcel;

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

my @maxgrad;
my @rmsgrad;
my @zpe;

FN2:
while ($fn = <$filelist>) {
    chomp $fn;
    my $fileout="$fn.out";
    open my $R, '<', $fileout or goto FN2;
    print "parsing $fileout\n";
    push @systemname, $fn;

    my @tmpgrad;
    my $tmpen;
    my $tmprms;
    my $tmpmax;
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
            my @tmp = split(/\s+/);
            $tmpen = $tmp[4];
            $tmpmax = $tmp[7];
            $tmprms = $tmp[9];
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
        if (/HARTREE\/MOLECULE/) {
		$_ =~ s/^\s+//;
		my ($tmp1,$tmp2,$tmp3,$tmp4) = split(/\s+/);
		push @zpe,$tmp1;
	}
            
    }

    push @energy, $tmpen;
    push @maxgrad, $tmpmax;
    push @rmsgrad, $tmprms;

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
print "Stationary Energy\t";
print "ZPE\t";
print "Maximum gradient\t";
print "RMS gradient\n";
my $jj=0;
foreach (@systemname) {
    print $systemname[$jj], "\t";
    print $energy[$jj], "\t";
    print $zpe[$jj], "\t";
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
print $W "Stationary Energy,";
print $W "ZPE,";
print $W "Maximum gradient,";
print $W "RMS gradient\n";
$jj=0;
foreach (@systemname) {
    print $W $systemname[$jj], ",";
    print $W $energy[$jj], ",";
    print $W $zpe[$jj], ",";
    print $W $maxgrad[$jj], ",";
    print $W $rmsgrad[$jj];
    print $W "\n";
    $jj++;
}

close $W;

#my $workbook = Excel:Writer::XLSX->new("$list_file_name.xlsx");
#my $worksheet = $workbook->add_worksheet();

#$worksheet->write(1,0,

