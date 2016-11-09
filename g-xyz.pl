#!/usr/bin/perl
#
# extract standard orientation
#

@sym = qw(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti Y Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn);

while(<>) {
    if (/^ NSERCH/) {$searchline = $_;}
    if (/COORDINATES OF ALL ATOMS ARE/) {
        print if $dbg;
        $_ = <>; $_ = <>;
        $i = 0;
        while(<>) {
            last if(/^\n/);
            $_ =~ s/^ //;
            chomp;
            $line[$i++] = $_;
        }
    }
}

print scalar(@line),"\n";
print $searchline;
for($k=0; $k < $i; $k++) {
    ($tmp1,$tmp2,$tmp3,$tmp4,$tmp5) = split /\s+/, $line[$k];
    print "$tmp1 $tmp3 $tmp4 $tmp5\n";
}
