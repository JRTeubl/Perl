#!/usr/bin/perl -w
use strict;

#/ifs/data/proteomics/tcga/samples/breast

my $dir="";
my $tcga ="";
my $uncid = "";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];}

open (OUT, ">>$dir/tcga_to_uncid.txt");

opendir (DIR, $dir);
my @allfolders = readdir DIR;
close DIR;
foreach my $folder (@allfolders){
    if ($folder=~/^(TCGA.+)/){
        $tcga = $1;
        print "$tcga\n";
        opendir (IN, "$dir/$folder");
        my @allfolder1 = readdir IN;
        close IN;
        foreach my $folder1 (@allfolder1){
            if ($folder1=~/rna/){
                opendir(RNA, "$dir/$folder/$folder1");
                my @allfolder2 = readdir RNA;
                close RNA;
                foreach my $folder2 (@allfolder2){
                    if ($folder2=~/fastq/){
                        opendir(FASTQ, "$dir/$folder/$folder1/$folder2");
                        my @allfolder3 = readdir FASTQ;
                        close FASTQ;
                        foreach my $folder3 (@allfolder3){
                            if ($folder3=~/^(UNCID_\d+)\./){
                                $uncid = $1;
                                print OUT "$uncid\t$tcga\n";
                            }
                        }
                    }
                }
            }
        }
    }else{print "Could not find TCGA folders\n"};
}