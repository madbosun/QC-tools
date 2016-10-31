#!/usr/bin/perl -w
  use strict;

  my $pid;
  my $option = '';
  my @shmids;
  my %userhash;
  my $user;

  $pid = open USR, "whoami |";
  $user = <USR>;
  $user =~ s/\s+//;
  close USR;

  if (@ARGV) {
      $option = $ARGV[0]
      if (scalar(@ARGV) > 1) or die "Only 1 input argument\n";
  }

  $pid = open IPCS, "ipcs -s |";

  while (<IPCS>) {
      my @dummy;
      chomp;
      s/^\s_//;
      if ($_ =~ /^key/) {
          #print "$_\n";
      }
      if ($_ =~ /^0/) {
          #print "$_\n";
          @dummy = split(/\s+/);
          push @shmids, $dummy[1];
          $userhash{$dummy[1]} = [] unless exists $userhash{$dummy[1]};
          $userhash{$dummy[1]} = $dummy[2];
      }
  }
  close IPCS;


if ($option =~ /dry/) {
    exit 0;
} else {
  print "Killing the following semaphore allocations:\n";
  foreach (@shmids) {
      if ($userhash{$_} =~ $user) {
          print "Killing $_\n";
          $pid = system("ipcrm -s $_");
      }
  }
  exit 0;
}
