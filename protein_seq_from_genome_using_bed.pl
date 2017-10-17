#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $filename="";
my $filename_vcf_germline="";
my $filename_vcf_somatic="";
my $dir=".";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $dir=$ARGV[1];} else { $dir="."; }
if ($ARGV[2]=~/\w/) { $filename_vcf_germline=$ARGV[2];} else { $filename_vcf_germline=""; }
if ($ARGV[3]=~/\w/) { $filename_vcf_somatic=$ARGV[3];} else { $filename_vcf_somatic=""; }


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
				
if ($error==0)
{
	my $filename_=$filename;
	$filename_=~s/\.bed$//;
	open (OUT,">$filename_.fasta");
	open (OUT_ALTPEP,">$filename.alt.pep.fasta");
	open (OUT_MOD,">$filename-mod.fasta");
	open (OUT_MODPEP,">$filename.mod.pep.fasta");
	open (BED,">$filename-mod.bed");
	open (LOG,">$filename-mod.log");
	open (STAT,">$filename-mod.stat");
	my $proteins_count=0;
	my $proteins_modified=0;
	my $count_stop_removed=0;
	my $count_stop_introduced=0;
	my $count_stop_removed_somatic=0;
	my $count_stop_introduced_somatic=0;
	my $protein_modifications_count=0;
	my @protein_modifications_distr=();
	my $protein_modifications_count_somatic=0;
	my @protein_modifications_distr_somatic=();
	my $count_variant_in_exon=0;
	my $count_variant_in_exon_somatic=0;
	my $count_variant_in_exon_old_error=0;
	my $count_variant_in_exon_old_error_somatic=0;
	my $count_variant_in_exon_nonsyn=0;
	my $count_variant_in_exon_nonsyn_somatic=0;
	my $line="";
	my %chr=();
	my %bed=();
	my %seq=();
	if (open (IN,"$filename"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^(\d+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)\s+(.+)$/)
			{
				my $chr=$1;
				#print "$chr\n";
				my $name=$5;
				#print "$name\n";
				$bed{$name}=$line;
				$chr{$chr}.="#$name#";
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	my %vcf_old=();
	my %vcf_new=();
	my %vcf_type=();
	my $germline_vcf=0;
	my $germline_only=0;
	my $somatic_only=0;
	my $both_vcf=0;
	my $both_vcf_differ=0;
	if (open (IN,"$filename_vcf_germline"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				my $chr="chr$1";
				my $pos=$2;
				my $id=$3;
				my $old=$4;
				my $new=$5;
				my $qaul=$6;
				$pos--;
				$new=~s/\,.*$//; #####
				$vcf_old{"$chr#$pos"}=$old;
				$vcf_new{"$chr#$pos"}=$new;
				$vcf_type{"$chr#$pos"}="G";
				$germline_vcf++;
				#print qq!$chr#$pos#$new\n!;
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	if (open (IN,"$filename_vcf_somatic"))
	{
		open (OUT_ONLY,">$filename_vcf_somatic-somatic_only.vcf");
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				my $chr="$1";
				my $pos=$2;
				my $id=$3;
				my $old=$4;
				my $new=$5;
				my $qaul=$6;
				$pos--;
				$new=~s/\,.*$//; #####
				if ($vcf_old{"$chr#$pos"}=~/\w/)
				{
					if ($vcf_new{"$chr#$pos"}!~/^$new$/)
					{
						#print LOG qq!$chr $pos: germline:$vcf_old{"$chr#$pos"}->$vcf_new{"$chr#$pos"} somatic:$old->$new\n!;
						$vcf_old{"$chr#$pos"}=$old;
						$vcf_new{"$chr#$pos"}=$new;
						$vcf_type{"$chr#$pos"}="S";
						print OUT_ONLY qq!$line\n!;
						$somatic_only++;
						$both_vcf_differ++;
					}
					$both_vcf++;
				}
				else
				{
					$vcf_old{"$chr#$pos"}=$old;
					$vcf_new{"$chr#$pos"}=$new;
					$vcf_type{"$chr#$pos"}="S";
					print OUT_ONLY qq!$line\n!;
					$somatic_only++;
				}
				#print qq!$chr#$pos#$new\n!;
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
		close(OUT_ONLY);
	}
	$germline_only=$germline_vcf-$both_vcf_differ;
	print qq!Somatic: $somatic_only\nGermline: $germline_only\nIn both VCF: $both_vcf (differ: $both_vcf_differ)\n!;
	
	foreach my $chr (sort keys %chr)
	{
		print qq!$chr\n!;
		if (open (IN,"$dir/chr.$chr.fa"))
		{
			print qq!opened $chr\n!;
			my $sequence="";
			$line=<IN>;
			chomp($line);
			my $chr_=$chr;
			$chr_=~s/^chr//i;
			if ($line=~/^>$chr\s*/ or $line=~/^>$chr_\s*/ or $line=~/\d+\s*/)
			{
				print"in line\n";
				while ($line=<IN>)
				{
					chomp($line);
					if ($line=~/^>/)
					{
						print qq!Error: > not expected: $line\n!;
					}
					else
					{
						$line=~s/\s+//g;
						if ($line!~/^[atcgATCGnN]+$/)
						{
							print qq!Error: unexpected character: $line\n!;
						}
						else
						{
							$sequence .= "\U$line";
						}
					}
				}
				my $temp=$chr{$chr};
				#print "$temp\n";
				while ($temp=~s/^#([^#]+)#//)
				{
					my $alt_splice=0;
					my $alt_splice1=0;
					my $alt_splice2=0;
					my $name=$1;
					print "$name\n\n";
					if ($name=~/^[^\-]+\-[^\-]+\-A\-([0-9]+)\-([0-9]+)/)
					{
						$alt_splice=1;
						$alt_splice1=$1;
						$alt_splice2=$2;
					}
					my $modified=0;
					print LOG qq!\n$name: $bed{$name}\n!;
					if ($bed{$name}=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
					{
						my $start=$2;
						#print "$start\n";
						my $end=$3;
						my $strand=$6;
						my $num=$10;
						my $segment_lengths="$11,";
						my $segment_starts="$12,";
						if ($strand=~/\+/) { if ($segment_lengths=~s/([0-9]+)\,$//) { my $temp=$1; $temp+=100; $segment_lengths.="$temp,"; } }
						else 
						{
							if ($segment_lengths=~s/^([0-9]+)\,//) { my $temp=$1; $temp+=100; $segment_lengths="$temp,$segment_lengths"; }
							if ($segment_starts=~s/^([0-9]+)\,//) { my $temp=$1; $temp-=100; $segment_starts="$temp,$segment_starts"; }
						}
						my $seq="";
						my $seq_original="";
						my $segment_starts_=$segment_starts;
						my $segment_lengths_=$segment_lengths;
						my %variants=();
						my %variants_=();
						my %variants2=();
						my %variants2_=();
						my $variant_loc=""; 
						my $seqment_lengths_splice=0;
						my $seqment_lengths_sum=0;
						my $seqment_count=0;
						my $segment_start_first=0;
						if ($segment_starts_=~/^([0-9\-]+)\,/) { $segment_start_first=$1; }
						while ($segment_starts_=~s/^([0-9\-]+)\,//)
						{
							my $segment_start=$1;
							if ($segment_lengths_=~s/^([0-9]+)\,//)
							{
								my $segment_length=$1;
								if ($alt_splice!=0 and $seqment_count<=$alt_splice1)
								{
									$seqment_lengths_splice+=$segment_length;
								}
								$seqment_lengths_sum+=$segment_length;
								#print qq!$name $seqment_count. $seqment_lengths_sum\n!;
								my $seq_=substr $sequence,$start+$segment_start,$segment_length;
								$seq_original.=$seq_; 
								print LOG qq!$start+$segment_start,$segment_length: $seq_\n!;
								for(my $i=0,my $j=$start+$segment_start;$i<$segment_length;$i++,$j++)
								{
									my $n=substr $seq_,$i,1;
									#print qq!$chr#$j: $n $vcf_new{"$chr#$j"}\n!;
									my $temp="";
									if ($vcf_new{"$chr#$j"}=~/\w/ and $vcf_new{"$chr#$j"}!~/^$vcf_old{"$chr#$j"}$/)
									{
										if ($n!~/^$vcf_old{"$chr#$j"}$/) 
										{ 
											$temp=qq! (Warning: $n\!\=$vcf_old{"$chr#$j"}) - Ignored!; 
											my $i_=$i+$segment_start_first+length($seq);
											print LOG qq!Variant $chr $j: $n,$vcf_old{"$chr#$j"}$temp->$vcf_new{"$chr#$j"} $i_\n!;
											$count_variant_in_exon_old_error++;
											if ($vcf_type{"$chr#$j"}=~/S/) 
											{
												$count_variant_in_exon_old_error_somatic++;
											}
										}
										else
										{
											$modified++;
											substr $seq_,$i,1,$vcf_new{"$chr#$j"};
											my $i_=$i+$segment_start_first+length($seq);
											print LOG qq!Variant $chr $j: $vcf_type{"$chr#$j"} $n,$vcf_old{"$chr#$j"}->$vcf_new{"$chr#$j"} $i_\n!;
											$count_variant_in_exon++;
											if ($vcf_type{"$chr#$j"}=~/S/) 
											{ 
												$count_variant_in_exon_somatic++; 
											}
											$variants{$i_}=qq!$vcf_old{"$chr#$j"}#$vcf_new{"$chr#$j"}#$vcf_type{"$chr#$j"}!;
											$variants2{$i_}=qq!$chr-$j!; 
										}
									}
								}
								$seq.=$seq_; 
								$seqment_count++;
							} else { print qq!Error parsing $bed{$name}\n!; }
						}
						if ($strand=~/\-/)
						{
							$seqment_lengths_splice=$seqment_lengths_sum-$seqment_lengths_splice;
						}
						my $name_=$name; $name_=~s/\-[^\-]+$//;
						#my $gene=$name; $gene=~s/^[^\-]+\-//;
						my $length=length($seq);
						my $protein="";
						my $protein_original="";
						my $description_="";
						if ($strand=~/\-/)
						{
							my $seq_ = reverse $seq;
							$seq=$seq_;
							$seq=~tr/ATCG/TAGC/;
							$seq_ = reverse $seq_original;
							$seq_original=$seq_;
							$seq_original=~tr/ATCG/TAGC/;
							foreach my $key (keys %variants) 
							{ 
								my $key_=$length+$segment_start_first-$key-1; 
								$variants_{$key_}=$variants{$key};
								$variants2_{$key_}=$variants2{$key};
								#print qq!$key->$key_\n!;
							}
						}
						else
						{
							foreach my $key (keys %variants)
							{
								$variants_{$key}=$variants{$key};
								$variants2_{$key}=$variants2{$key};
							}
						}
						$modified=0;
						my $modified_somatic=0;
						my $stop_found=0;
						my $triplet_count=0;
						for(my $n=0;$n<$length and $stop_found==0;$n=$n+3)
						{
							my $n_=$n+2;
							my $triplet = substr($seq_original, $n, 3);
							if (length($triplet)==3)
							{
								if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
								$protein_original.=$mapping{$triplet}; 
							}
							$triplet = substr($seq, $n, 3);
							if (length($triplet)==3)
							{
								if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
								$protein.=$mapping{$triplet}; 
								my $triplet_old=$triplet;
								my $triplet_type="";
								for (my $i=$n, my $j=0; $i<=$n_; $i++,$j++) 
								{
									if ($variants_{$i}=~/^([^#]+)#([^#]+)#([^#]+)$/)
									{
										my $old=$1;
										$triplet_type=$3;
										if ($strand=~/\-/) { $old=~tr/ATCG/TAGC/; } 
										substr $triplet_old,$j,1,$old;
										$variant_loc=$variants2_{$i};
									}
								}
								if ($mapping{$triplet_old}!~/[\w\*]/) { $mapping{$triplet_old}="X"; }
								#print qq!$name $triplet_type $n-$n_: $triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet}\n!;
								if ($mapping{$triplet}=~/\*/)
								{
									$stop_found=1;
								}
								if ($mapping{$triplet}=~/\*/ and $mapping{$triplet_old}!~/\*/)
								{
									print LOG      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet}\n!;
									#$description_.=qq!$variant_loc-$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet},!;
									$description_.=qq!$variant_loc-$triplet_type $n-$n_-$triplet_old->$triplet - $mapping{$triplet_old}$triplet_count$mapping{$triplet},!;
									$modified++;
									$count_variant_in_exon_nonsyn++;
									$count_stop_introduced++;
									if ($triplet_type=~/S/) 
									{
										$modified_somatic++;
										$count_variant_in_exon_nonsyn_somatic++;
										$count_stop_introduced_somatic++;
									}
								}
								else
								{
									if ($mapping{$triplet}!~/\*/ and $mapping{$triplet_old}=~/\*/)
									{
										print LOG      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet}\n!;
										#$description_.=qq!$variant_loc-$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet},!;
										$description_.=qq!$variant_loc-$triplet_type $n-$n_-$triplet_old->$triplet - $mapping{$triplet_old}$triplet_count$mapping{$triplet},!;
										$modified++;
										$count_variant_in_exon_nonsyn++;
										$count_stop_removed++;
										if ($triplet_type=~/S/) 
										{
											$modified_somatic++;
											$count_variant_in_exon_nonsyn_somatic++;
											$count_stop_removed_somatic++;
										}
									}
									else
									{
										if ($mapping{$triplet}!~/^$mapping{$triplet_old}$/)
										{
											print LOG      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet}\n!;
											#$description_.=qq!$variant_loc-$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count:$mapping{$triplet_old}->$mapping{$triplet},!;
											$description_.=qq!$variant_loc-$triplet_type $n-$n_-$triplet_old->$triplet - $mapping{$triplet_old}$triplet_count$mapping{$triplet},!;
											$modified++;
											$count_variant_in_exon_nonsyn++;
											if ($triplet_type=~/S/) 
											{
												$modified_somatic++;
												$count_variant_in_exon_nonsyn_somatic++;
											}
										}
									}
								}
							}
							$triplet_count++;
						}
						$proteins_count++;
						$protein_modifications_count+=$modified;
						$protein_modifications_distr[$modified]++;
						$protein_modifications_count_somatic+=$modified_somatic;
						$protein_modifications_distr_somatic[$modified_somatic]++;
						$protein=~s/\*$//;
						$protein_original=~s/\*$//;
						if (length($protein_original)>6)
						{
							$protein_original=~s/^([^\*]+)\*.*$/$1/;
							print OUT qq!>$name \n$protein_original\n!;
							if ($alt_splice!=0)
							{
								my $k=int($seqment_lengths_splice/3);
								my $temp1=substr $protein_original,0,$k;
								$temp1=~s/^.*\*([^\*]*)$/$1/;
								my $cleavage=0;
								my $temp1_="";
								while($temp1=~s/(.)(.)$//)
								{
									my $aa1=$1;
									my $aa2=$2;
									if ($cleavage<2)
									{
										$temp1_="$aa2$temp1_";
										$temp1.="$aa1";
										if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
										{
											$cleavage++;
										}
									}
								}
								if ($cleavage==0) { $temp1_=""; }
								my $count_RK = $temp1=~tr/RK/RK/;
								if ($temp1!~s/^[^RK]*[RK]//) { $temp1=""; }
								for(my $l=0;$l<$count_RK-2;$l++) { $temp1=~s/^[^RK]*[RK]/$1/; }
								my $temp2=substr $protein_original,$k;
								$temp2=~s/([^\*]*)\*.*$/$1/;
								$cleavage=0;
								my $temp2_="";
								while($temp2=~s/^(.)(.)//)
								{
									my $aa1=$1;
									my $aa2=$2;
									if ($cleavage<2)
									{
										$temp2_.="$aa1";
										$temp2="$aa2$temp2";
										if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
										{
											$cleavage++;
										}
									}
								}
								if ($cleavage<2) { $temp2_.=$temp2; }
								if ($cleavage==0) { $temp2_=""; }
								if (length($temp1_)>0 and length($temp2_)>0 and length($temp1_)+length($temp2_)>6)
								{
									print OUT_ALTPEP qq!>$name-bridge \n$temp1_$temp2_\n!;
								}
							}
						}
						if ($modified!=0)
						{
							$proteins_modified++;
							$protein=~s/^([^\*]+)\*.*$/$1/;
							my $temp=$description_;
							while($temp=~s/^([^\,]+),\s*//)
							{
								my $temp1=$1; 
								if ($temp1=~s/\- ([A-Z]+)([0-9]+)([A-Z]+)$//) #added s to take out the values
								{
									my $aa_old=$1;
									my $loc=$2;
									my $aa=$3;
									#print "$aa_old$loc$aa\t$temp1\n";
									print LOG qq!---$1:$2-$3\n!;
									my $temp1_= substr $protein,0,$loc;
									my $aa_= substr $protein,$loc,1;
									my $temp2_= substr $protein,$loc;
									if ($temp1_=~/[RK]$/) { $temp1_=""; }
									$temp1_=~s/^.*[RK]([^P])/$1/;
									$temp2_=~s/([RK])[^P].*$/$1/;
									print LOG qq!--$temp1_\n--$temp2_\n!;
									if ($temp2_=~s/^$aa//) 
									{
										if ($aa=~/\*/)
										{
											if (length($temp1_)>6)
											{
												#print OUT_MODPEP qq!>$name-$temp1_\_$aa_old-stop $temp1- $temp1_($aa_old->$aa)\n$temp1_\n!;
												print OUT_MODPEP qq!>$name $temp1 $aa_old$loc$aa-stop\n$temp1_\n!;
											}
											#$description_.=" $temp1_($aa_old->$aa)";
											#$description_.=" $aa_old$loc$aa"; 
										}
										else
										{
											if (length("$temp1_$aa$temp2_")>6)
											{
												#print OUT_MODPEP qq!>$name-$temp1_\_$aa_old-$aa\_$temp2_ $temp1- $temp1_($aa_old->$aa)$temp2_\n$temp1_$aa$temp2_\n!;
												print OUT_MODPEP qq!>$name $temp1 $aa_old$loc$aa\n$temp1_$aa$temp2_\n!;
											}
											#$description_.=" $temp1_($aa_old->$aa)$temp2_";
											#$description_.=" $aa_old$loc$aa"; 
										}
									}
									else
									{
										$description_.=" $temp1_-$aa_old->Error($aa)-$temp2_";
									}
								}
								#$description_.=",";
							}
							
							if (length($protein)>6)
							{
								print LOG qq!>$name-variant $description_\n$protein\n!;
								print OUT_MOD qq!>$name-variant $description_\n$protein\n!;
								
								my $protein_length=length($protein);
								my $segment_starts_=$segment_starts;
								my $segment_lengths_=$segment_lengths;
								my $segment_starts__="";
								my $segment_lengths__="";
								my $segment_count=0;
								my $length_sum=0;
								my $length_sum_=0;
								my $end_found=0;
								my $start_=$start;
								my $end_=$end;
								my $segment_start_first=0;
								if ($strand=~/\+/)
								{
									while ($segment_starts_=~s/^([0-9\-]+)\,//)
									{
										my $segment_start=$1;
										if ($segment_lengths_=~s/^([0-9]+)\,//)
										{
											my $segment_length=$1;
											$length_sum+=$segment_length;
											if ($end_found==0)
											{
												if ($length_sum>3*$protein_length)
												{
													$end_found=1;
													$segment_length-=$length_sum-3*$protein_length;
													$segment_length+=3;
												}
												$segment_starts__.="$segment_start,";
												$segment_lengths__.="$segment_length,";
												$end_=$start+$segment_start+$segment_length;
												$segment_count++;
											}
										}
									}
								}
								else
								{
									my $segment_start_first=0;
									if ($segment_starts_=~/^([0-9\-]+)\,/) { $segment_start_first=$1; }
									while ($segment_starts_=~s/([0-9\-]+)\,$//)
									{
										my $segment_start=$1;
										if ($segment_lengths_=~s/([0-9]+)\,$//)
										{
											my $segment_length=$1;
											$length_sum+=$segment_length;
											if ($end_found==0)
											{
												my $segment_end=$segment_start+$segment_length;
												if ($length_sum>3*$protein_length)
												{
													$end_found=1;
													$segment_length-=$length_sum-3*$protein_length;
													$segment_length+=3;
												}
												$segment_start-=$segment_start_first;
												$segment_starts__="$segment_start,$segment_starts__";
												$segment_lengths__="$segment_length,$segment_lengths__";
												$start_=$start+$segment_end-$segment_length;
												$segment_count++;
											}
										}
									}
								}
								$segment_lengths__=~s/\,$//;
								$segment_starts__=~s/\,$//;
								print BED qq!$chr\t$start_\t$end_\t$name-variant\t1000\t$strand\t$start_\t$end_\t0\t$segment_count\t$segment_lengths__\t$segment_starts__\n!;
								print LOG qq!---$name\t$chr\t$start_\t$end_\t$name-variant\t1000\t$strand\t$start_\t$end_\t0\t$segment_count\t$segment_lengths__\t$segment_starts__\n!;
							}
						}
					} else { print qq!Error parsing $bed{$name}\n!; }
				}
			} else { print qq!Error in name $chr: $line\n!; }
			
			close(IN);
		}
	}
	print STAT qq!
		$proteins_count proteins
		$proteins_modified proteins modified
		$count_stop_removed stop codons removed
		$count_stop_introduced stop codons introduced
		Somatic: $count_stop_removed_somatic stop codons removed
		Somatic: $count_stop_introduced_somatic stop codons introduced
		$protein_modifications_count total modifications
		$protein_modifications_count_somatic somatic modifications
		Somatic: $somatic_only
		Germline: $germline_only
		In both VCF: $both_vcf (differ: $both_vcf_differ)
		$count_variant_in_exon variants in exons
		$count_variant_in_exon_old_error variants in exons where old does not match genome
		$count_variant_in_exon_nonsyn non-synonymous variants
		Somatic: $count_variant_in_exon_somatic variants in exons
		Somatic: $count_variant_in_exon_old_error_somatic variants in exons where old does not match genome
		Somatic: $count_variant_in_exon_nonsyn_somatic non-synonymous variants
	!;
	print STAT qq!\nnumber of modifications\tnumber of proteins\tnumber of proteins (somatic variants)\n!;
	#for (my $k=0;$k<=100;$k++) { print STAT qq!$k\t$protein_modifications_distr[$k]\t$protein_modifications_distr_somatic[$k]\n!; }
	close(OUT);
	close(LOG);
	close(STAT);
	close(BED);
}
