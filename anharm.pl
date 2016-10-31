#!/usr/bin/perl

#
# extract Gaussian Anharmonic data
#

my $go=0;

my @AnSpectra;
my @AnIntensity;
my @HarmSpectra;
my @HarmIntensity;

while(<>) {

    if (/Anharmonic Infrared Spectroscopy/) { $go = 1; }
    if ($go > 0) {
	if (/Fundamental Bands/) { $go = 2; }
	if (/Overtones/) { $go = 3; }
	if (/Combination Bands/) { $go = 4; }

	if ($go > 1) {
	    $_ = <>;
	    $_ = <>;
	    while(<>) {
		if (/^\n/) {
		    $go = 1;
		    last
		}
		$_ =~ s/^\s+//;
		chomp;
		my @tmp = split(/\s+/);

                if ($go == 2) {
		    push @AnSpectra, $tmp[2];
		    push @AnIntensity, $tmp[4];
		    push @HarmSpectra, $tmp[2];
		    push @HarmIntensity, $tmp[4];
		} elsif ($go == 3) {
		    push @AnSpectra, $tmp[2];
		    push @AnIntensity, $tmp[3];
		} elsif ($go == 4) {
		    push @AnSpectra, $tmp[3];
		    push @AnIntensity, $tmp[4];
		}

		#print "$_\n";
	    }
	}
    }

}

#print " ### CUT HERE FOR _AN_HARMONIC SPECTRA ####\n";
print "#System\n#Method\n";

for($k=0; $k < scalar(@AnSpectra); $k++) {
    printf "%-10.3f %14.7f\n", $AnSpectra[$k], $AnIntensity[$k];
}

#print "\n ### CUT HERE FOR HARMONIC SPECTRA ####\n";
#print "#System\n#Method\n";
#
#for($k=0; $k < scalar(@HarmSpectra); $k++) {
#    printf "%-10.3f %14.7f\n", $HarmSpectra[$k], $HarmIntensity[$k];
#}
