#!/usr/bin/perl -w
use strict;
#################################################################
#
# Author: Jason Byrd
# Contact: madbosun@gmail.com
# date: July 24 2006
# version 1.2
#
#################################################################
#
# Read in a pdb file, and take a spherical chunk out of the system
#   based on $Asearch.  The radius of the sphere is defined by 
#   $rcutoff, an input variable set by the flag -r. 
#################################################################
#
## define variables here
# work with solution as global res
  my $solres = "WAT";
  my $molres = "MOL";
  my $rcutoff;
# we are working with H2O here, so... and we count O's here
  my $slide = 2;  my $Asearch = "O";
#help subroutine
  sub help {
    print "\nPDB ion system reduction script\n";
    print "\nsyntax: pickclean [options] <filename>\n";
    print "\noptions:\n";
    print "  -s res  : solute residue name\n";
    print "  -w res  : solvent residue name\n";
    print "  -r #    : radial cutoff\n";
    print "\n";
    exit 1;
  }
#
#################################################################

  use Getopt::Long;
  GetOptions("s:s" => \$molres,"w:s" => \$solres, "r:s" => \$rcutoff);
  my $infilename=$ARGV[0];

  if (!$rcutoff) {
    print "WARNING: no radial cutoff option set.\n";
    help();
  } elsif (@$ARGV[0]) {
    print "WARNING: no input file name.\n";
    help();
  }
  print "Using radial cutoff $rcutoff\n";
  if ($Asearch =~ "O") {
    print "Working with H20 solvate, with residue name $solres, change line 24 if other solvate is wanted\n";
  }

##open the file
  $infilename =~ s/.pdb//;
  open(F,"$infilename.pdb") || print "bad file input name\n",help();
  my @line=<F>;
  close(F);
  open(W,">picked-$rcutoff-$infilename.com");
  open(Z,">picked-$rcutoff-$infilename.zmat");
  print W "#p test\n";
  print W " \n";
  print W "picked $infilename with $rcutoff radius\n";
  print W " \n";
  print W "0,1\n";
  print Z "picked $infilename with $rcutoff radius\n";

## start reading the pdb file in, and parse line by line
  my $molcount=0;

# parse out central monomer 
# center on first atom in that monomer
  my @ligx; my @ligy; my @ligz;
  print "Centering on atom 1 of residue $molres\n";
  my $cx=0; my $cy=0; my $cz=0;
    my $jj=0;
    foreach (@line) {
      s/^\s+//;
      if ($_ =~ /^ATOM/) {
        $jj++;
        my ($t1,$t2,$t3,$t4,$t5,$t6,$t7,$t8,$t9,$t10) = split(/\s+/);
        chomp;
	$t3 =~ s/\d+$//;
	$t3 =~ s/W$//;
          if ($t4 =~ $molres) {
	    $molcount++;
	    if ($jj == 1) {
		$cx = -1*$t6; $cy = -1*$t7; $cz = -1*$t8;
        	my $cx = $cx; my $cy = $cy; my $cz = $cz;
	    }
	    my $x=$t6+$cx; my $y=$t7+$cy; my $z=$t8+$cz;
	    push(@ligx,$x); push(@ligy,$y); push(@ligz,$z);
            print W sprintf ("%-4s",$t3);
            print W sprintf ("%13s",sprintf("%.*f",6,$x));
            print W sprintf ("%13s",sprintf("%.*f",6,$y));
            print W sprintf ("%13s",sprintf("%.*f",6,$z));
            print W "\n";
            print Z sprintf ("%-4s",$t3);
            print Z sprintf ("%13s",sprintf("%.*f",6,$x));
            print Z sprintf ("%13s",sprintf("%.*f",6,$y));
            print Z sprintf ("%13s",sprintf("%.*f",6,$z));
            print Z "\n";
          } elsif ($t4 =~ $solres) {
	      last;
          }
      }
    }

# center and store coordinates
  print "Centering global system\n";
  my $iswater=0;
  my $atomcount=$molcount;
  my $solcount=1;
  my $r;
    my $q;
    foreach (@line) {
      s/^\s+//;
      if ($_ =~ /^ATOM/) {
        chomp;
        my ($t1,$t2,$t3,$t4,$t5,$t6,$t7,$t8,$t9,$t10) = split(/\s+/);
	my $x=$t6+$cx;
	my $y=$t7+$cy;
	my $z=$t8+$cz;
	$t3 =~ s/\d+$//;
	$t3 =~ s/W$//;
	my $pstring=sprintf ("%-4s",$t3);
	$pstring .= sprintf ("%13s",sprintf("%.*f",6,$x));
	$pstring .= sprintf ("%13s",sprintf("%.*f",6,$y));
	$pstring .= sprintf ("%13s",sprintf("%.*f",6,$z));
	$r = 1000.0;
        for (my $ii = 0; $ii < scalar(@ligx); $ii++) {
	    $q = sqrt( ($x-$ligx[$ii])*($x-$ligx[$ii]) + ($y-$ligy[$ii])*($y-$ligy[$ii]) + ($z-$ligz[$ii])*($z-$ligz[$ii]));
	    if ($q < $r){$r = $q;}
	}
        if ($t4 =~ $solres) {
	  if ($t3 =~ $Asearch && $r <= $rcutoff && $iswater == 0) {
            print W $pstring;
            print W "\n";
            print Z $pstring;
            print Z "\n";
	    $iswater=1;
	    $atomcount++;
	    $solcount++;
	  } elsif ($iswater == 1) {
            print W $pstring;
            print W "\n";
            print Z $pstring;
            print Z "\n";
	    $atomcount++;
	    $iswater=2;
	  } elsif ($iswater == 2) {
            print W $pstring;
            print W "\n";
            print Z $pstring;
            print Z "\n";
	    $atomcount++;
	    $iswater=0;
	  }
	}
    }
    }

    print "Now there are $atomcount atoms with $solcount monomers\n";

    print Z "\n";
    print Z "*ACES2(calc=ccsd\n";
    print Z "basis=FOOBASIS\n";
    print Z "cc_conv=7\n";
    print Z "scf_expstart=3,scf_exporder=20,scf_conv=6\n";
    print Z "estate_sym=1,estate_tol=4\n";
    print Z "spherical=on,coordinates=cartesian,symmetry=none)\n";
    print Z "\n";
    print Z "*SIP\n";
    print Z "MAXMEM=120000\n";
    print Z "SIAL_PROGRAM = MOI_rhf.siox\n";
    print Z "SIAL_PROGRAM = scf_frag_lowmem.siox\n";
    print Z "SIAL_PROGRAM = mcpt2_corr.siox\n";
    print Z "\n";
    print Z "*FRAGLIST\n";
    print Z "$solcount R_DISP R_ELST\n";
    print Z "$molcount ";

    print "assuming water solvent: 3 atoms per solvent monomer\n";
    for (my $ii = 0; $ii < $solcount - 1; $ii++) {
	print Z "3 ";
    }
    print Z "\n";
    for (my $ii = 1; $ii <= $molcount; $ii++) {
	print Z "$ii ";
    }
    print Z "\n";

    for (my $ii = $molcount + 1; $ii <= $atomcount; $ii += 3) {
	print Z $ii," ",$ii+1," ",$ii+2,"\n";
    }
    print Z "\n";

    print W "\n";

    close(W);
    close(Z);
