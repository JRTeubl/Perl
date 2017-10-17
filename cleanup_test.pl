
my $dir="";
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }

if ($dir=~/test/i)
{
	system(qq!rm $dir/*.mgf!);
	system(qq!rm $dir/*.xml!);
	system(qq!rm $dir/*.sh!);
}