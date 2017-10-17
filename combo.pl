#!/usr/bin/perl -w
use strict;

my %pep = ();
my $pep="";
my $dir="";
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }

open (OUT, ">>$dir/combo.txt");

if (opendir(DIR,"$dir")){
        my @allfiles=readdir DIR;
        closedir DIR;
        foreach my $filename (@allfiles){
            if ($filename=~/\.out$/i) {
                print qq!$filename\n!;
		open (IN,"$dir/$filename") || die "Could not open $dir/$filename\n";
                while (my $line = <IN>){
                    if ($line =~ /^([A-Z]+)\s+(.+)$/){
                        $pep = $1;
                        #print "$pep\n";
                        $pep{$pep}= 1;
                    }
                }
            }
        }
} else {die "could not open $dir\n"};

print OUT "$_\n" for keys %pep;
