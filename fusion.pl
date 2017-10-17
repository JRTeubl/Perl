#!/usr/bin/perl -w
# 

#maps between ensemble and refseq protein identifiers using human_hgnc_refseq_ensembl.txt file
#this file contains symbol \t name \t refseqID \t ensemblID \n
#also use the hgnc_ENSB_NP_Accession.txt to convert accession numbers to ensemble genes and refseq IDs
#and Homo_sapiens_ensembl_list.txt to convert between ensembl genes, transcripts and proteins 
use strict;
no warnings; 

my $file = "";
my @allfiles=();
my $error = 0;
my %fusion=();
my $name="";
my $sample=0;
my $chr1="";
my $chr2="";
my $gene1="";
my $gene2="";
my $loc1="";
my $loc2="";
my $seq="";
my $root=""; 

if ($ARGV[0]=~/\w/) { $file=$ARGV[0];} else {$error=1;}
if ($ARGV[1]=~/\w/) { $root=$ARGV[1];} else {$error=1;}

my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
				"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
				"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
				"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
				
				"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
				"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
				"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
				"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
				
				"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
				"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
				"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
				"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
				
				"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
				"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
				"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
				"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");

if ($error == 0)
{
    if (open(IN, "$file"))
    {
	while (my $line =<IN>)
	{
	    if ($line=~/<TABLE CELLPADDING\=3 BORDER\=\"1\"><TR><TD ALIGN\=\"LEFT\"><H1>UNCID\_([0-9]+)</)
	    {
		$name=$1;
		$name="UNCID_$name";
		$sample=1;
		$seq=""; 
	    }
	    if ($line=~/<TD ALIGN\=\"LEFT\"><H1>([A-Z]*[0-9]*)<\/H1><\/TD>/ && $sample==2)
	    {
		$gene2=$1; 
	    }
	    if ($line=~/<TD ALIGN\=\"LEFT\"><H1>chr([0-9]+)<\/H1><\/TD>/ && $sample==2)
	    {
		$chr2=$1;
		$chr2="chr$chr2";
	    }
	    if ($line=~/<TD ALIGN\=\"RIGHT\"><H1>([0-9]+)<\/H1><\/TD>/ && $sample==2)
	    {
		$loc2=$1;
		$sample=3; 
	    }
	    if ($line=~/<TD ALIGN\=\"LEFT\"><H1>([A-Z]*[0-9]*)<\/H1><\/TD>/ && $sample==1)
	    {
		$gene1=$1; 
	    }
	    if ($line=~/<TD ALIGN\=\"LEFT\"><H1>chr([0-9]+)<\/H1><\/TD>/ && $sample==1)
	    {
		$chr1=$1;
		$chr1="chr$chr1";
	    }
	    if ($line=~/<TD ALIGN\=\"RIGHT\"><H1>([0-9]+)<\/H1><\/TD>/ && $sample==1)
	    {
		$loc1=$1;
		$sample=2; 
	    }
	    if ($line=~/<H3>left flanking sequence \- ([ATCG]+)<\/H3>/ && $sample==3)
	    {
		$seq=$1; 
	    }
	    if ($line=~/<H3>right flanking sequence \- ([ATCG]+)<\/H3>/ && $sample==3)
	    {
		$seq.=$1;
		$fusion{$name}.="#$gene1-$chr1-$loc1-$gene2-$chr2-$loc2#$seq#";
		$sample=0; 
	    }
	}
    }
    else
    {
	print "Unable to open fusion result.hml file\n"; 
    }
    close IN;

    open (OUT, ">$root/fusion.fasta") or die "could not open output file!";
    open (OUT2, ">$root/proteome.fasta") or die "can't open file";
    {
	foreach my $key (keys %fusion)
	{
	    my $info=$fusion{$key};
	    while($info=~s/^#([^#]*)#([^#]*)#//)
	    {
		my $label=$1;
		my $sequence=$2;
		$name=">$label-$key"; 
		#translate the sequence
		if ($name=~/\w/ && $sequence=~/\w/)
		{
		    foreach my $direction ("fwd","rev")
		    {
			my $seq="";
			my $size=length($sequence);
			if ($direction =~ /^fwd$/) { $seq=$sequence; } 
			else 
			{ 
			    $seq = reverse $sequence;
			    $seq =~ tr/ATCG/TAGC/;
			}
			for (my $k=0;$k<3;$k++)
			{#for each reading frame:
			    my $protein="";
			    for(my $n=$k;$n<$size;$n=$n+3)
			    {
				my $triplet = substr($seq, $n, 3);
				if ($mapping{$triplet} =~ /[\w\*]/) { $protein .= $mapping{$triplet}; } # '*' is stop codon
				else { $protein.="X"; } # X is unknown, doesn't code for anything must be error in sequence
			    }
			    my $temp=$protein;
			    my $longest_orf="";
			    if ($temp=~/\*/)
			    {
				while ($temp =~ s/^([^\*]*)\*//) 
				{
				    my $orf_seq=$1;
				    if (length($orf_seq)>6)
				    {
					print OUT2 qq!>$name\_$direction\_fr$k $direction frame $k\n$orf_seq\n!;
				    }
				}
			    }
			    else
			    {
				print OUT2 qq!>$name\_$direction\_fr$k $direction frame $k\n$temp\n!;
			    }
			    print OUT qq!>$name\_$direction\_fr$k $direction frame $k\n$protein\n!;
			}
		    }
		}
		
	    }
	    
	}
    }
}
else
{
    print "ERROR!"; 
}

