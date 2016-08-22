#!/usr/bin/perl
#
# extract GAMESS hess block
#

while(<>) {
    if (/\$HESS/) {
	$line[0] = $_;
        $i = 1;
        while(<>) {
            $line[$i++] = $_;
            last if(/\$END/);
        }
    }
}

for($k=0; $k < $i; $k++) {
print $line[$k]
}
