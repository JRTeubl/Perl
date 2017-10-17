#!/usr/bin/perl -w
use strict;

my $tcga = "";
if ($ARGV[0]=~/\w/) {$tcga=$ARGV[0]};
open (OUT, ">>new_cccd.txt");
open (TCGA, $tcga);
while (my $line = <TCGA>){
    $line=~s/\n/\t/gs;
    if($line=~s/(^\t\t){25}/$1\n/gs){
        print OUT $line;
    }
}


