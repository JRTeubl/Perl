#!/usr/bin/perl -w
use strict;

my $tcga = "";
my $list1 = "";
my $list2 = "";
my %id = ();
my $id = "";


if ($ARGV[0]=~/\w/) {$tcga=$ARGV[0]};
if ($ARGV[1]=~/\w/) {$list1=$ARGV[1]};
if ($ARGV[2]=~/\w/) {$list2=$ARGV[2]};
open (OUT, ">>new_cccd.txt");


#foreach (keys %id){
#    print "$_\n";
#}


open (TCGA, $tcga);
while (my $line = <TCGA>){
    if ($line=~/^(\w+\-.{2}-.{4})\t(.*)$/){
    $id = $1;
    #print "$id\n";
    $id{$id}= $2;
    }else{print OUT "$line"};
    
}

open (IN1, $list1);
while (my $line = <IN1>){
    if ($line=~/^(\w+\-.{2}-.{4}).*$/){
        $id = $1;
        print OUT "$id\t$id{$id}\n";
        #print "$1\n";
        #if (defined $id{$id}){
        #    print OUT "$id\t$id{$id}\n";
        #}
        #print "$id{$1}\n";
    }
}

open (IN2, $list2);
while (my $line = <IN2>){
    if ($line=~/(.+)/){
        $id=$1;
        #print "$1\n";
        if (exists $id{$id}){
            print OUT "$id\t$id{$id}\n";
        }
    }
}
#for my $key (keys %id){
#    if ($id{$key}!= 1){
#        print OUT "$key\t$id{$key}\n"
#    }
#}


#while (my ($k, $v)= each %id){
#    print "$k\t$v\n";
#}

#if (exists $id{$id}){
#           $id{$id}= $2;
#            print "hash exists";
#        }
