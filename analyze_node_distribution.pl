my %stat=();
my %stat_error=();
my %unique_models=();
my %unique_models_count=();
my $count=0;
my $count_=0;
my $count_done=0;

my $dir="";
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }

if (opendir(dir,"$dir"))
{
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $filename (@allfiles)
	{
		if ($filename=~/^(.*)\.o([0-9]+)$/i)
		{
			my $base=$1;
			my $number=$2;
			if (open(IN,"$dir/$filename"))
			{
				my $node="";
				while($line=<IN>)
				{
					if ($line=~/^Running on\s+(\w+)\s*$/)
					{
						$node=$1;
						$stat{$node}++;
						if (open(IN_,"$base.e$number"))
						{
							my $line_=<IN_>;
							if ($line_=~/\w/) { $stat_error{$node}++; }
							close(IN_);
						}
						$count++;
					}
					if ($line=~/^Unique models = ([0-9]+)/)
					{
						$unique_models{"$node"}+=$1;
						$unique_models_count{"$node"}++;
						$count_done++;
					}
				}
				close(IN);
				$count_++;
			}
		}
	}
	print qq!\n$dir\n$count_ o files\n$count "Running on..."\n$count_done done\n!;
	
	print qq!node\tjobs\tdone\tavg\terror\n!;
	foreach my $node (sort keys %stat)
	{
		if ($unique_models_count{"$node"}>0) { $unique_models{"$node"}/=$unique_models_count{"$node"}; }
		print qq!$node\t$stat{$node}\t$unique_models_count{"$node"}\t$unique_models{"$node"}!;
		if ($stat_error{$node}>0)
		{
			print qq!\t***$stat_error{$node} errors***\n!;
		}
		else
		{
			print qq!\t\n!;
		}
			
	}
}