#!/usr/bin/perl -w
use strict;

my $file="";
my %pep=();
if ($ARGV[0]=~/\w/) { $file=$ARGV[0];} else { print "no file found\n";}
open (IN, "$file");
open (OUT, ">>rat_ensm_pep_clean.txt");

while (my $line = <IN>){
    if ($line =~ /^(.+)$/){
        $pep{$1}=1;
    }
}

foreach my $keys (keys %pep){
    print OUT "$keys\n";
}
