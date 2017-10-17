#!/usr/bin/perl
#
use strict;

my $error=0;
my $dir="";
my $expect_threshold="";
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }
if ($ARGV[1]=~/\w/) { $expect_threshold=$ARGV[1];} else { $expect_threshold=1e-2; }

if ($error==0)
{
	$dir=~s/\\/\//g;
	open (LOG,qq!>$dir.$expect_threshold.log!) || die "Could not open output\n";
	my $line="";
	my $filename ="";
	if (opendir(DIR,"$dir"))
	{	
		my @allfiles=readdir DIR;
		closedir DIR;
		my $output= "peptide_list.out";
			
		foreach $filename (@allfiles)
		{	open (IN,"$dir/$filename") || die "Could not open $dir/$filename\n";
			open (OUT,qq!>>$output!) || die "Could not open output\n";
			if ($filename=~/\.t\.xml$/i)
			{       print qq!$filename\n!;	
				my $reversed=1;
				my $pep="";
				my $expect="";
				my $peptide_proteins="";
				my $name = "";
				while ($line=<IN>)
				{	
					if ($line=~/^\<protein\s+.*label="([^\"]+)"/)
					{
						my $protein_name=$1;
						my $protein=$protein_name;
						$protein=~s/^(\S+)\s.*$/$1/;
						if ($protein_name!~/\:reversed$/) { $reversed=0; }
						$peptide_proteins.="#$protein#";
					}
					if ($line=~/^\<domain\s+id="([0-9\.edED\+\-]+)".*expect="([0-9\.edED\+\-]+)".*mh="([0-9\.edED\+\-]+)".*delta="([0-9\.edED\+\-]+)".*seq="([A-Z]+)"/)
					{
						my $expect_=$2;
						my $pep_=$5;
						
						$pep_=~tr/L/I/;
						if ($expect!~/\w/ or $expect_<=$expect) { $expect=$expect_; $pep=$pep_; }
						#print "$pep\n";
					}
					if($line=~/<note label=\"description\">(.+)\spep.+$/ or $line=~/<note label=\"description\">(.+)\|(.+)\|.+$/){
						$name = $2;
						#print "$name\n";
					}
					if($line=~/<note label=\"Description\">(.+?)<\/note>/)	
					{
						my $title=$1;
						#print "$title\n";
						chomp($title);
						$title=~s/\n//;
						$title=~s/\r//;
						if ($reversed==0 and $expect<$expect_threshold)
						{
							print OUT qq!$pep\t$name\n!;
						}
						$reversed=1;
						$pep="";
						$expect="";
						$peptide_proteins="";
					}
				}
					
				close(IN);
				close(OUT);	
			}
		}
	}
	close(LOG);
}