#!/usr/bin/perl
# USAGE:
# unix command line:
use strict;
use warnings;

my $inFile = shift @ARGV;    # read file names from terminal input
open IN, '<', $inFile or die "Cannot open '$inFile' because: $!";



while (<IN>){
chomp;
next if ($_ =~ m/^##/);
if ($_ =~ /^#/){
print "$_\n";
}
else {
my @tmpo=split("\t", $_);
my $count=0;
for (my $i=9; $i<=$#tmpo; $i++){
my @data=split("\:", $tmpo[$i]);
next if ($data[2] eq ".");
if ($data[2] < 10) {
$count++;
}#if
}#for
next if ($count > 0);
print "$_\n";
}# else
}

close IN;



