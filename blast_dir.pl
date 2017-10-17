#!/usr/local/bin/perl

use strict;

my $error=0;
my $dir="";
if ($ARGV[0]=~/\w/) { $dir="$ARGV[0]";} else { $dir="."; }

if ($error==0)
{
	system(qq!formatdb -p F -i genome_rat.fasta!);
	if (opendir(DIR,"$dir"))
	{
		my @allfiles=readdir DIR;
		closedir DIR;
		foreach my $filename (@allfiles)
		{
			if ($filename=~/\.fa$/i)
			{
				print qq!$filename\n!;
				system(qq!formatdb -p F -i $filename!);
				system(qq!blastall -p blastn -d genome_rat.fasta -i $filename -o /ifs/data/proteomics/projects/rat/1600_iTRAQ/Output/blast/$filename.out!);
			}
		}
	}
}

