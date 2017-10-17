#!/usr/bin/perl -w
use strict;

my %kinase = ();
my @data = ();
my $prot = "";
open(KIN, "/Users/jrteubl/Desktop/kinase_list.txt");
while(my $line = <KIN>){
    chomp $line;
    print "$line\n";
    $kinase{$line}= 1;
}
close KIN;

open(DATA, "/Users/jrteubl/Desktop/phospho_gene_dataframe.txt");
while(my $line = <DATA>){
    if ($line=~/^(\d+|\w+)\t.+$/){
        $prot = $1;
        #print "$prot\n";
        push @data, $prot if exists $kinase{$prot};
    }
}

open(OUT, ">/Users/jrteubl/Desktop/kinases_to_use.txt");
foreach my $item (@data){
    print OUT "$item\n";
}