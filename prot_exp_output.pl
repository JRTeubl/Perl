#!/usr/bin/perl -w
use strict;


my $dir = "";
my $clinical = "";
my $id = "";
my $subtype = "";
my %store = ();
my $barcode = "";
my $name = "";
my $exp = "";

if ($ARGV[0]=~/\w/) {$dir=$ARGV[0]};
if ($ARGV[1]=~/\w/) {$clinical = $ARGV[1]};

my $newdir= "$dir/output";
mkdir $newdir;

open (CLIN, $clinical) || die "Could not open $clinical\n";
while (my $line = <CLIN>){
    if ($line =~/((.+)\t*){30}/){ ###fix this line
        $id = $1;
        print "$id\n";
        $subtype = $21;
        $store{$id}= $subtype;
    }
}

opendir (DIR, "$dir");
my @allfiles = readdir DIR;
close DIR;
foreach my $file (@allfiles){
    if ($file =~/^(.+)__(.+)__(.+)\-(.+)\-(.+)\-(.+)__(.+)$/){
        $barcode = "$3-$4-$5";
        open (OUT, ">>$newdir/$barcode.output.txt");
        open (IN, "$dir/$file");
        while (my $line = <IN>){
            if ($line =~/(.+)\t(.+)\t(.+)\t(.+)/){
                $name = $3;
                $exp = $4;
                if ($name ne "gene name"){
                    foreach my $key (keys %store){
                        print $key;
                        if ($key eq $barcode){
                            print OUT "$key\t$store{$key}\n";
                            print OUT "$name\t$exp\n"
                        }
                    }
                }
            }
        }
    }
}