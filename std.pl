#!/usr/bin/perl
#
# extract standard orientation
#

@sym = qw(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti Y Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn);

while(<>) {
    if (/Standard orientation:/) {
        print if $dbg;
        $_ = <>; $_ = <>; $_ = <>; $_ = <>;
        $i = 0;
        while(<>) {
            last if(/-------------------/);
            ($ian,$x[$i],$y[$i],$z[$i]) = (split)[1,3,4,5];
            $s[$i++] = $sym[$ian];
            printf "%-4s %12s %12s %12s\n", $sym[$ian], $x, $y, $z if $dbg;
        }
    }
}

print "struct\n";
print "$i\n";
for($k=0; $k < $i; $k++) {
    printf "  %-4s %12s %12s %12s\n", $s[$k], $x[$k], $y[$k], $z[$k];
}
print "\n";
