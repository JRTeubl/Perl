
my $dir=".";
my $species="human";
my $precursor_mass_error_ppm=20;
my $fragment_mass_error_ppm=60;
my $fixed_mods="";
my $potential_mods="";
my $param_dir="/ifs/data/proteomics/programs/tandem/bin";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }
if ($ARGV[1]=~/\w/) { $species=$ARGV[1];} else { $species="human"; }
if ($ARGV[2]=~/\w/) { $precursor_mass_error_ppm=$ARGV[2];} else { $precursor_mass_error_ppm=20; }
if ($ARGV[3]=~/\w/) { $fragment_mass_error_ppm=$ARGV[3];} else { $fragment_mass_error_ppm=60; }
if ($ARGV[4]=~/\w/) { $fixed_mods=$ARGV[4];} else { $fixed_mods="57.022@C"; }
if ($ARGV[5]=~/\w/) { $potential_mods=$ARGV[5];} else { $potential_mods=""; }
if ($ARGV[6]=~/\w/) { $param_dir=$ARGV[6];} else { $param_dir="/ifs/data/proteomics/programs/tandem/bin"; }

my $result_dir="$species-$precursor_mass_error_ppm-$fragment_mass_error_ppm-";
$result_dir.=GetDateTime_();
mkdir "$dir/$result_dir";
open(OUT_SH,">$dir/$result_dir/tandem.sh");

print OUT_SH qq!#\!/bin/bash\n#\$ -cwd\ncd $dir/$result_dir\n!;

if (opendir(dir,"$dir"))
{
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $filename (@allfiles)
	{
		if ($filename=~/^(.*)\.mgf$/i or $filename=~/^(.*)\.mzML$/i)
		{
			my $filename_=$1;
			open(OUT,">$dir/$result_dir/$filename_\_input.xml");
			print OUT qq!<?xml version="1.0"?>
<bioml>
	<note type="input" label="list path, default parameters">$param_dir/default_input.xml</note>
	<note type="input" label="list path, taxonomy information">$param_dir/taxonomy.xml</note>
	<note type="input" label="spectrum, fragment monoisotopic mass error">$fragment_mass_error_ppm</note>
	<note type="input" label="spectrum, parent monoisotopic mass error plus">$precursor_mass_error_ppm</note>
	<note type="input" label="spectrum, parent monoisotopic mass error minus">$precursor_mass_error_ppm</note>
	<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
	<note type="input" label="spectrum, fragment monoisotopic mass error units">ppm</note>
	<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
	<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
!;
if ($fixed_mods=~/\w/)
{
	print OUT qq!		<note type="input" label="residue, modification mass">$fixed_mods</note>\n!;
}

if ($potential_mods=~/\w/)
{
	print OUT qq!	<note type="input" label="residue, potential modification mass">$potential_mods</note>\n!;
}
print OUT qq!
	<note type="input" label="protein, taxon">$species</note>
	<note type="input" label="spectrum, path">../$filename</note>
	<note type="input" label="output, path">$filename_.xml</note>
</bioml>
!;
			close(OUT);
	
			print OUT_SH qq!qsub ./$filename_.sh\n!;
	
			open(OUT,">$dir/$result_dir/$filename_.sh");
			print OUT qq!#\!/bin/bash
#\$ -cwd
echo "Running on `hostname`" 
module load tandem
tandem.exe ./$filename_\_input.xml
!;
			close(OUT);
		}
	}
}

close(OUT_SH);
system(qq!chmod u+x $dir/$result_dir/*.sh!);

system(qq!$dir/$result_dir/tandem.sh!);

sub GetDateTime_
{
	my $sec="";
	my $min="";
	my $hour="";
	my $mday="";
	my $mon="";
	my $year="";

	($sec,$min,$hour,$mday,$mon,$year) = localtime();

	if ($sec<10) { $sec="0$sec"; }
	if ($min<10) { $min="0$min"; }
	if ($hour<10) { $hour="0$hour"; }
	if ($mday<10) { $mday="0$mday"; }
	$mon++;
	if ($mon<10) { $mon="0$mon"; }
	$year+=1900;
	$date="$year-$mon-$mday-$hour-$min-$sec";
	
	return $date;
}
