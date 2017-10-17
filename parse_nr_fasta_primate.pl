#!/usr/bin/perl -w
use strict;


my $error=0;
my $dir = "/ifs/data/proteomics/projects/Siegel/Species_Id/databases/current/nr_20151110";
my $fasta = "/ifs/data/proteomics/projects/Siegel/Species_Id/databases/current/nr_20151110/nr.fasta";
my $file_name = "";
if ($ARGV[0]=~/\w/) { $fasta = $ARGV[0]; }
#if ($ARGV[1]=~/\w/) { $dir = $ARGV[1]; } else { $dir = "."; }
if ($fasta=~/(\w+)(\.fasta)$/) { $file_name = $1; }

my $gi_num = "";
my %name = ();
my $names = "";
my $seq = "";

my %names=();
if (open(IN, "/ifs/data/proteomics/projects/Siegel/Species_Id/databases/current/nr_20151110/taxonomy/names.dmp"))
{
	while(my $line = <IN>)
	{
		if ($line=~/^\s*([0-9]+)\s+\|\s+(.*)$/)
		{
			my $node=$1;
			my $name=$2;
			$name=~s/\t.\t.*$//;
			#print qq!#$node#$name#\n!;
			$names{$name}=$node;
		}
	}
	close(IN);
}
my %tree=();
if (open(IN, "/ifs/data/proteomics/projects/Siegel/Species_Id/databases/current/nr_20151110/taxonomy/nodes.dmp"))
{
	while(my $line = <IN>)
	{
		if ($line=~/^\s*([0-9]+)\s+\|\s+([0-9]+)/)
		{
			#print qq!$1 $2\n!;
			$tree{$1}=$2;
		}
	}
	close(IN);
}
my %primate=();
foreach my $node (keys %tree)
{
	my $done=0;
	if (exists $primate{$node})
	{
		if ($primate{$node}==1)
		{
			$done=1;
		}
	}
	if ($done==0)
	{
		$primate{$node}=0;
		my $parent=$node;
		while($parent!=1)
		{
			$parent=$tree{$parent};
			if ($parent==9443)
			{
				$primate{$node}=1;
			}
		}
	}
}
if (open(OUT, ">$dir/primates.txt"))
{
	print OUT qq!name\tid\tprimates\n!;
	foreach my $name (sort keys %names)
	{
		print OUT qq!$name\t$names{$name}\t$primate{$names{$name}}\n!;
	}
	close(OUT);
}

if (open(IN, "$fasta"))
{
	if (open(OUT, ">$dir/$file_name-edited.fasta"))
	{
		open(OUT1, ">$dir/species_table.txt");
		open(LOG, ">$dir/parse_nr_fasta.log");
		open(OUT_NO, ">$dir/$file_name-edited-not-primate.fasta");
		open(OUT1_NO, ">$dir/species_table-not-primate.txt");
		open(LOG_NO, ">$dir/parse_nr_fasta-not-primate.log");
		my $gi_num_previous="";
		my $line_previous="";
		while (my $line = <IN>)
		{
			chomp($line);
			#print "$line\n";
			if ($line =~ /^>gi\|(\d+)\|/)
			{
				$gi_num = "$1";
				#print "$gi_num\n";
				if ($line_previous=~/\w/)
				{
                                        $line_previous =~ s/' '//g;
					my @temp = split(/\b/, $line_previous);
					my $primate=0;
					%name = ();
                                        my $item_prev = "";
                                        my %spec_genus = ();
                                        my @name_list = ();
					foreach my $item (@temp)
					{
						if ($item =~/^\s*$/){
                                                    next
                                                }else{
                                                    if ($item_prev=~/^\s*\[/){
                                                    $spec_genus{$item}= 1;
                                                    }
                                                    if (exists($spec_genus{$item_prev})){
                                                        push @name_list, "$item_prev $item";
                                                    }
                                                    $item_prev = $item;
                                                }
                                        }
					my @name_list_vert = ();
                                        foreach my $item_name (@name_list){
						print "$item_name\n";
						if (exists $names{$item_name})
						{
							if (exists $primate{$names{$item_name}})
							{
								if ($primate{$names{$item_name}}==1)
								{
									push @name_list_vert, $item_name;
									$primate=1;
								}
							}
						}  
						else 
						{ 
						    if ($primate==1) { print LOG qq!Error: $gi_num_previous $item_name\n!; } else { print LOG_NO qq!Error: $gi_num_previous $item_name\n!; } 
						}
                                        }
                                        my @name_uniq = uniq(@name_list_vert);
					$names = join(',', @name_uniq);
					if ($primate==1)
					{
						print OUT1 "$gi_num_previous\t$names\n";
						if ($gi_num_previous=~/\w/)
						{
							if ($seq!~/[^A-Za-z]/)
							{
								print OUT ">$gi_num_previous\n$seq\n";
							}
							else
							{
								print LOG "Sequence error: $gi_num_previous: $seq\n";
							}
						}
					}
					#else
					#{
					#	print OUT1_NO "$gi_num_previous\t$names\n";
					#	if ($gi_num_previous=~/\w/)
					#	{
					#		if ($seq!~/[^A-Za-z]/)
					#		{
					#			print OUT_NO ">$gi_num_previous\n$seq\n";
					#		}
					#		else
					#		{
					#			print LOG_NO "Sequence error: $gi_num_previous: $seq\n";
					#		}
					#	}
					#}
				}
				$gi_num_previous=$gi_num;
				$line_previous=$line;
				$seq="";
			}
			else
			{
				$seq.=$line
			}
		}
        }
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

