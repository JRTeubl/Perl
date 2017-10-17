#!/usr/bin/perl -w
use strict;

#this program separates each transcript sequence to it's own file so each file can be blasted against the rat genome

#file = sequences.txt

my $file="";
my $dir="";
if ($ARGV[0]=~/\w/) {$file="$ARGV[0]";}else{ print "cannot find file\n";}
if ($ARGV[1]=~/\W/) {$dir = "$ARGV[1]";}

open (IN, "$file");

my $line;
my $ensm;
my $seq;
my $header;
while ($line = <IN>){
    if ($line=~/^(>)(ENS\w{4}\d{11})\s+(.+)$/){ 
        $header = $1;
        $ensm = $2;
        #print "$ensm\n";
        $seq = $3;
        open (OUT,">$dir/$ensm.fa");
        print OUT "$header$ensm\n$seq\n";
    }else{print "no match\n";}
}
