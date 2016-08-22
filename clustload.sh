#!/usr/bin/perl -w

  use strict;
  use warnings qw(void syntax severe redefine);
  use warnings FATAL => qw(io overflow);
  no warnings qw(uninitialized);

  my @stuff;
  my $ps;
  my $pid;
  my $memory=0.0;
  my $cpu=0.0;
  my $submemory=0.0;
  my $subcpu=0.0;
  my $search;
  my $count=0;
  my @shmids;
  my %userhash;
  my $user;


# check for input, take only first argument regardless of any thing
  if (@ARGV) {
          $search = $ARGV[0];
          if (scalar(@ARGV) > 1) {
                  print "Only 1 input argument\n";
                  exit 1;
          }
  }

# general pid for ps pipe
  $pid = open PS, "ps -o user,pmem,nlwp,pcpu,cputime,comm,pid -e |";

# general loop over the input pipe
  LINE: while(<PS>) {
    chomp;
    s/^\s+//;
    @stuff = split(/\s+/);
# print words on the first line then skip analysis
        if ($count == 0) {
                $count++;
            if ($stuff[0] =~ /^USER/i) {&psprint (@stuff);}
                next LINE;
    }
# we are counting total memory here
    $memory = $memory + $stuff[1];
    $cpu = $cpu + $stuff[3];
        if ($search) {
                next LINE if $stuff[5] !~ $search;
                &psprint (@stuff);
                $submemory = $submemory + $stuff[1];
                $subcpu = $subcpu + $stuff[3];
        }
# search for 10 percent cpu
        if ($stuff[3] > 10.0) {&psprint (@stuff);}
# search for 5 percent memory
        elsif ($stuff[1] > 5.0) {&psprint (@stuff);}
  }

  sub psprint {
      print sprintf("%-10s",$stuff[0]);
      print sprintf("%-6s",$stuff[1]);
      print sprintf("%-5s",$stuff[2]);
      print sprintf("%-5s",$stuff[3]);
      print sprintf("%-11s",$stuff[4]);
      print sprintf("%-18s",$stuff[5]);
      print sprintf("%-6s",$stuff[6]);
      print "\n";
  }
  close PS;
  if ($search) {
          print "\nSpecified Program Memory Usage: $submemory%\n";
          print "Specified Program CPU Usage: $subcpu%\n";
          }

#  $pid = open IOSTAT, "iostat |";
#
#  print "\n---For Hex---\n";
#  while (<IOSTAT>) {
#         chomp;
#         if ($_ =~ /^avg-cpu/) {
#             print "$_\n";
#         } elsif ($_ =~ /^         /) {
#             print "$_\n";
#         }
#         if ($_ =~ /^Device:/) {
#             print "\n$_\n";
#         } elsif ($_ =~ /^sda1/) {
#             print "$_\n";
#         } elsif ($_ =~ /^sdc1/) {
#             print "$_\n";
#         }
#
#  }
#  print " \n";
#  close(IOSTAT);

  $pid = open IPCS, "ipcs -m |";

  while (<IPCS>) {
      my @dummy;
      chomp;
      s/^\s_//;
          print "$_\n";
      if ($_ =~ /^key/) {
      }
      if ($_ =~ /^0/) {
      @dummy = split(/\s+/);
      push @shmids, $dummy[1];
      $userhash{$dummy[1]} = [] unless exists $userhash{$dummy[1]};
      $userhash{$dummy[1]} = $dummy[2];
      }
  }

  close IPCS;

  $pid = open FREE, "free -m -t |";
  while (<FREE>) {
          chomp;
          print "$_\n";
  }
  close FREE;
  print "\nTotal ps memory use: $memory%\n";
  print "Total cpu use: $cpu%\n";
  exit 0;

  $pid = open PROCMEM, "cat /proc/meminfo |";
  my $printer;
  while (<PROCMEM>) {
          chomp;
          if ($_ =~ /^MemTotal:/){print "$_\n";}
          if ($_ =~ /^MemFree:/){print "$_\n";}
          if ($_ =~ /^Active:/){print "$_\n";}
          if ($_ =~ /^Inactive:/){print "$_\n";}
  }
  close PROCMEM;
  print "\nTotal ps memory use: $memory%\n";
  print "Total cpu use: $cpu%\n";
  exit 0;
