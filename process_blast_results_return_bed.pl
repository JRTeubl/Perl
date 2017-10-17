#!c:/perl/bin/perl.exe
# Parses blast output to modify the sequences with least number of gaps, insertions and deletions.
# This program replaces the sequences with best matches. 

use strict;
my $error=0;
my $dir="";
if ($ARGV[0]=~/\w/) { $dir="$ARGV[0]";} else { $dir="."; }

if ($error==0)
{
	if (opendir(DIR,"$dir"))
	{
		my @allfiles=readdir DIR;
		closedir DIR;
		foreach my $filename (@allfiles)
		{
			if ($filename=~/\.out$/i)
			{
				print qq!$filename\n!;
				open (OUT,">>parse_blast.bed");
				if (open (IN,"$filename"))
				{
					my $start=0;
					my $end=0;
					my $seq="";
					my %start_query=();
					my %seq_query=();
					my %seq_query_=();
					my %end_query=();
					my $start_=0;
					my $end_=0;
					my $seq_="";
					my %start_sbjct=();
					my %seq_sbjct=();
					my %seq_sbjct_=();
					my %end_sbjct=();
					my $line="";
					my $line_="";
					my $count=0;
					my $chr="";
					my $chr_="";
					my %chr=();
					my $strand="";
					my $trim_num=3;
					my $percentage_cutoff=80;
					my @to_sort_chr=();
					my %results=();
					my %segment_count=();
					my %segment_start_q=();
					my %segment_seq_q=();
					my %segment_end_q=();
					my %segment_start_s=();
					my %segment_seq_s=();
					my %segment_end_s=();
					
					
					while ($line=<IN>)
					{
						chomp($line);
						#print qq!$line\n!;
						if ($line=~/^>(\d+\s+dna).+$/) 
						{
							$chr=$1;
							$chr_=$1;
							$chr{$chr}=1;
							#print qq!$chr\n! 
						}
						if ($line=~/^\s+Strand\s=\s(.*)/)
						{
							$strand=$1;
							$chr="$chr_-$strand";
							#print "$strand\n";
						}
						else
						{ 
							if ($line=~/^Query:\s*([0-9]+)\s+(\S*)\s+([0-9]+)/)
							{
								$start=$1;
								$seq=$2;
								$end=$3;
								$start_query{"$chr#$count"}=$start;
								$end_query{"$chr#$count"}=$end;
								$seq_query_{"$chr#$start#$end"}=$seq; 
								$line=<IN>;
								$line=<IN>;
								if ($line=~/^Sbjct:\s*([0-9]+)\s+(\S*)\s+([0-9]+)/)
								{
									$start_=$1;
									$seq_=$2;
									$end_=$3;
									$start_sbjct{"$chr#$count"}=$start_;
									$end_sbjct{"$chr#$count"}=$end_;
									$seq_sbjct_{"$chr#$start_#$end_"}=$seq_; 
									my $k=0;
									
									
									
									my $k_max=length($seq_query_{"$chr#$start#$end"});
									my $match_percentage=0;   
									while($k<$k_max)
									{
										my $amino_q = uc(substr($seq_query_{"$chr#$start#$end"}, $k, 1));
										my $amino_s = uc(substr($seq_sbjct_{"$chr#$start_#$end_"}, $k, 1)); 
										if($amino_q eq $amino_s) { $match_percentage++; }
										$k++; 
									} 
									if(($match_percentage/$k_max)*100 >= $percentage_cutoff) 
									{ 
										$k=0;
										my $ok=0;
										while($ok==0)
										{
											my $comp_q = uc(substr($seq_query_{"$chr#$start#$end"}, $k, $trim_num)); 
											my $comp_s = uc(substr($seq_sbjct_{"$chr#$start_#$end_"}, $k, $trim_num)); 
											if ($comp_q eq $comp_s) { $ok=1; } else { $k++; }
										} 
										if($ok==1)
										{
											$start_query{"$chr#$count"}+=$k;
											$start_sbjct{"$chr#$count"}+=$k; 
											$seq_query_{qq!$chr#$start_query{"$chr#$count"}#$end!}=substr $seq_query_{"$chr#$start#$end"},$k; 
											$seq_sbjct_{qq!$chr#$start_sbjct{"$chr#$count"}#$end_!}=substr $seq_sbjct_{"$chr#$start_#$end_"},$k; 
										}
										$k=-$trim_num;
										$ok=0;
										while($ok==0)
										{
											my $comp_q = uc(substr($seq_query_{qq!$chr#$start_query{"$chr#$count"}#$end!},$k,$trim_num)); 
											my $comp_s = uc(substr($seq_sbjct_{qq!$chr#$start_sbjct{"$chr#$count"}#$end_!},$k,$trim_num)); 
											if ($comp_q eq $comp_s) { $ok=1; } else { $k--; }
										}
										if($ok==1)
										{
											$end_query{"$chr#$count"}+=$k+$trim_num;
											$end_sbjct{"$chr#$count"}+=$k+$trim_num; 
											if($k==-$trim_num) 
											{ 
												$seq_query{qq!$chr#$start_query{"$chr#$count"}#$end_query{"$chr#$count"}!}=substr $seq_query_{qq!$chr#$start_query{"$chr#$count"}#$end!},0;   
												$seq_sbjct{qq!$chr#$start_sbjct{"$chr#$count"}#$end_sbjct{"$chr#$count"}!}=substr $seq_sbjct_{qq!$chr#$start_sbjct{"$chr#$count"}#$end_!},0;  
											}
											else 
											{ 	
												$seq_query{qq!$chr#$start_query{"$chr#$count"}#$end_query{"$chr#$count"}!}=substr $seq_query_{qq!$chr#$start_query{"$chr#$count"}#$end!},0,$k+$trim_num; 
												$seq_sbjct{qq!$chr#$start_sbjct{"$chr#$count"}#$end_sbjct{"$chr#$count"}!}=substr $seq_sbjct_{qq!$chr#$start_sbjct{"$chr#$count"}#$end_!},0,$k+$trim_num;
											}
											if ($start_query{"$chr#$count"}+$trim_num-1<=$end_query{"$chr#$count"}) 
											{
												#print qq!$count. ($match_percentage out of $k_max)\n$seq_query{"$chr#$start_query{\"$chr#$count\"}#$end_query{\"$chr#$count\"}"} ($start_query{"$chr#$count"}-$end_query{"$chr#$count"})\n$seq_sbjct{"$chr#$start_sbjct{\"$chr#$count\"}#$end_sbjct{\"$chr#$count\"}"}\n!;
											}
										}
										if ($start_query{"$chr#$count"}+$trim_num-1<=$end_query{"$chr#$count"}) { $count++; }
									}
									else 
									{ 
										$seq_query_{"$chr#$start#$end"}=""; $seq_sbjct_{"$chr#$start_#$end_"}=""; #print qq!less than cutoff ($match_percentage out of $k_max)\n!; 
									}
								}
							}
						}
					}
					close(IN);
					
					foreach my $chr_ (sort keys %chr)
					{
						foreach my $strand ("Plus / Plus", "Plus / Minus")
						{
							$chr="$chr_-$strand";
							#print qq!\n\n\n$chr\n!;
							my %list=();
							my $query_sequence="";
							my $sbjct_sequence="";
							for(my $i=0;$i<$count;$i++)
							{
								$query_sequence=$seq_query{qq!$chr#$start_query{"$chr#$i"}#$end_query{"$chr#$i"}!};
								$sbjct_sequence=$seq_sbjct{qq!$chr#$start_sbjct{"$chr#$i"}#$end_sbjct{"$chr#$i"}!};
								$list{qq!$start_query{"$chr#$i"}#$query_sequence#$end_query{"$chr#$i"}#$start_sbjct{"$chr#$i"}#$sbjct_sequence#$end_sbjct{"$chr#$i"}!}=qq!$start_sbjct{"$chr#$i"}!;
							}
							my @sorted_list=();
							my $iterator=0;
							foreach my $value (sort {$list{$a} <=> $list{$b}} keys %list)
							{
								my $temp=$value; 
								if($temp=~/([^#]+)#([^#]+)#([^#]+)#([^#]+)#([^#]+)#([^#]+)/)
								{
									$sorted_list[$iterator]="$1#$2#$3#$4#$5#$6"; 
									$iterator++;
								}
							}
							my $segment_count=-1;
							my $match=0;
							my $match1=0;
							my $compare_s=0;
							my $compare_e=0;
							my $compare1_s=0;
							my $compare1_e=0;
							for(my $i=0;$i<@sorted_list;$i++)
							{   
								if($sorted_list[$i]=~/([^#]+)#([^#]+)#([^#]+)#([^#]+)#([^#]+)#([^#]+)/)
								{ 
									my $start_q=$1;
									my $seq_q=uc($2);
									my $end_q=$3;
									my $start_s=$4;
									my $seq_s=uc($5);
									my $end_s=$6;
									if ($segment_count==-1 or $segment_end_s{"$chr#$segment_count"} < $start_s -1)
									{
										$segment_count++;
										#print qq!\n\n!;
										$segment_start_q{"$chr#$segment_count"}=$start_q;
										$segment_seq_q{"$chr#$segment_count"}=$seq_q;
										$segment_end_q{"$chr#$segment_count"}=$end_q;
										$segment_start_s{"$chr#$segment_count"}=$start_s;
										$segment_seq_s{"$chr#$segment_count"}=$seq_s;
										$segment_end_s{"$chr#$segment_count"}=$end_s;	
									}
									else
									{   
										my $aminoacid_q="";
										my $aminoacid_s="";
										my $j=$start_s-$segment_start_s{"$chr#$segment_count"} + 1; 
										my $temp=substr($segment_seq_s{"$chr#$segment_count"},0,$j);
										my $insert_count=($temp=~tr/-/-/);
										$j+=$insert_count; 
										my $j_max=0;
										if ($segment_end_s{"$chr#$segment_count"}<=$end_s) { $j_max=length($segment_seq_s{"$chr#$segment_count"})+1; } else { $j_max=$j+length($seq_s); }
										#print qq!$j,$j_max\n!;
										$match=0;
										$compare_s=$j;
										$compare_e=$j_max;
										########### printing overlapped regions
										if($j != $j_max)
										{
											my $overlap_s1=substr $segment_seq_s{"$chr#$segment_count"},$j-1,$j_max; 
											my $overlap_q=substr $segment_seq_q{"$chr#$segment_count"},$j-1,$j_max;
											#print qq!q :$overlap_q\ns1:$overlap_s1\n!; 
										} 
										########### printing overlapped regions
										while($j<$j_max)
										{
											$aminoacid_q = uc(substr($segment_seq_q{"$chr#$segment_count"}, $j, 1));
											$aminoacid_s = uc(substr($segment_seq_s{"$chr#$segment_count"}, $j, 1)); 
											if($aminoacid_s eq $aminoacid_q) { $match++; }
											$j++;
										}
										if ($segment_end_q{"$chr#$segment_count"}<=$end_q) { $j_max=$segment_end_q{"$chr#$segment_count"} - $start_q +1; } else { $j_max=$end_q - $start_q +1; }
										$j=0; 
										$match1=0; 
										$compare1_s=$j;
										$compare1_e=$j_max;
										########### printing overlapped regions
										if($j != $j_max) 
										{ 
											my $overlap_s=substr $seq_s,$j,$j_max;
											#print qq!s2:$overlap_s\n!; 
										} 
										########### printing overlapped regions
										while($j<$j_max)
										{
											$aminoacid_q = uc(substr($seq_q, $j, 1));
											$aminoacid_s = uc(substr($seq_s, $j, 1));
											if($aminoacid_s eq $aminoacid_q) { $match1++; }
											$j++;
										}
										if($match > $match1)
										{ 
											if($segment_end_q{"$chr#$segment_count"}<=$end_q)
											{
												$segment_seq_q{"$chr#$segment_count"} .= substr $seq_q,$segment_end_q{"$chr#$segment_count"} - $start_q + $insert_count;
												$segment_seq_s{"$chr#$segment_count"} .= substr $seq_s,$segment_end_q{"$chr#$segment_count"} - $start_q + $insert_count;
												$segment_end_q{"$chr#$segment_count"} = $end_q;
												$segment_end_s{"$chr#$segment_count"} = $end_s; 
											}
											else
											{	
											} 
										} 
										else
										{ 
											if($segment_end_q{"$chr#$segment_count"}<=$end_q)
											{
												$segment_seq_q{"$chr#$segment_count"} = substr $segment_seq_q{"$chr#$segment_count"},0,$start_q - $segment_start_q{"$chr#$segment_count"} + $insert_count;
												$segment_seq_s{"$chr#$segment_count"} = substr $segment_seq_s{"$chr#$segment_count"},0,$start_q - $segment_start_q{"$chr#$segment_count"} + $insert_count; 
												$segment_seq_q{"$chr#$segment_count"} .= $seq_q;
												$segment_end_q{"$chr#$segment_count"} = $end_q;
												$segment_seq_s{"$chr#$segment_count"} .= $seq_s;
												$segment_end_s{"$chr#$segment_count"} = $end_s; 
											}
											else
											{
												my $temp = substr $segment_seq_q{"$chr#$segment_count"},0,$start_q - $segment_start_q{"$chr#$segment_count"} + $insert_count;
												$temp .= $seq_q;
												$temp .= substr $segment_seq_q{"$chr#$segment_count"},$end_q+1,$segment_end_q{"$chr#$segment_count"} - $end_q + $insert_count;
												$segment_seq_q{"$chr#$segment_count"}=$temp;
												$temp="";
												$temp = substr $segment_seq_s{"$chr#$segment_count"},0,$start_q - $segment_start_q{"$chr#$segment_count"} + $insert_count;
												$temp .= $seq_s;
												$temp .= substr $segment_seq_s{"$chr#$segment_count"},$end_q+1,$segment_end_q{"$chr#$segment_count"} - $end_q + $insert_count;
												$segment_seq_s{"$chr#$segment_count"}=$temp; 
											} 
										}
									}
									#print qq!$segment_count.$i ($match,$match1) ($compare_s-$compare_e,$compare1_s-$compare1_e).\n$segment_seq_q{"$chr#$segment_count"}\n$segment_seq_s{"$chr#$segment_count"}\n!;
								} 		
							}
							#print OUT qq!$chr\n!;
							my $matches=0;
							for(my $i=0;$i<=$segment_count;$i++)
							{
								$matches+=abs($segment_end_q{"$chr#$i"}-$segment_start_q{"$chr#$i"}+1);
								#print OUT qq!$i. $segment_start_q{"$chr#$i"}-$segment_end_q{"$chr#$i"}: $segment_seq_q{"$chr#$i"}\n$i. $segment_start_s{"$chr#$i"}-$segment_end_s{"$chr#$i"}: $segment_seq_s{"$chr#$i"}\n!;
								$results{"$chr#$i#start"}=qq!$segment_start_q{"$chr#$i"}!;
								$results{"$chr#$i#end"}=qq!$segment_end_q{"$chr#$i"}!;
							}
							@to_sort_chr=(@to_sort_chr,"$matches#$chr");
							$segment_count{$chr}=$segment_count;
							
						}
					}
					my @sorted_chr = sort {$b<=>$a} @to_sort_chr;
					$chr_=$sorted_chr[0];
					$chr_=~s/^([0-9]+)#//;
					if ($chr_=~/^(.*)\-(.*)$/)
					{
						my $chr= $1;
						#print "$chr\n";
						my $strand = $2;
						#print "$strand\n";
						if ($strand =~/Plus\s\/\sPlus/)
						{
							#print "$strand\n";
							my $g = 1;
                                                        my %length = ();
							my %offset = ();
                                                        my %exons = ();
							my %store= ();
							my %seq_length = ();
							my $i=0;
							my $r = $i+1;
							my $stop = 1;
							#print "segment count\t";
							#print $segment_count{$chr_};
							while ($i<=$segment_count{$chr_})
							{
								while ($stop == 1)
								{
									$store{qq!$segment_start_s{"$chr_#$i"}#$segment_end_s{"$chr_#$i"}#$segment_seq_s{"$chr_#$i"}!}= $g;
									my $l=$results{"$chr_#$i#end"}-$results{"$chr_#$i#start"}+1;
									my $o=$results{"$chr_#$i#start"}-$results{"$chr_#0#start"};
									my $len = length $segment_seq_q{"$chr_#$i"};
									$seq_length{$g} += $len;
									$length{$g}.="$l,";
									$offset{$g}.="$o,";
									$exons{$g} += 1;
									if ($segment_end_s{"$chr_#$i"} <= $segment_start_s{"$chr_#$r"} and $segment_end_q{"$chr_#$i"} <= $segment_start_q{"$chr_#$r"})
									{
										$i++;
										$r++;
									}else{
										$stop = 0;
									}
								}
								
								$i++;
								$r++;
                                                                $g++;
								$stop = 1;
							}
                                                        #print OUT "$i, $r, $g\n";
							my @store=();
							foreach my $key (sort {$a <=> $b} keys %store)
                                                        {
                                                            push @store, "$store{$key}\.$key\n";
                                                        }
							#print OUT @store;
                                                        
                                                        my @length =();
							foreach my $key (sort {$seq_length{$b} <=> $seq_length{$a}} keys %length)
                                                        {
                                                            push @length, "$key\t$seq_length{$key}\n";
                                                        }
                                                        #print OUT @length;
                                                        my $line = $length[0];
                                                        my $group ="";
                                                        my $length="";
                                                        if ($line =~/^(\d+)\t(\d+)/)
                                                        {
                                                            $group = $1;
							    #print "$group\n";
                                                            $length = $2;
                                                        }
                                                        my $newstart = "";
                                                        my $newend = "";
                                                        my $sequence = "";
                                                        my %final = ();
                                                        my $p ="";
                                                        for ($p = 0; $p<$i; $p++)
                                                        {
                                                            #print OUT "$p\n";
                                                            if ($store[$p]=~/^($group)\.(\d+)#(\d+)#(.+)$/)
                                                            {
                                                                $newstart = $2;
                                                                #print OUT "start = $newstart\n";
                                                                $newend = $3;
                                                                #print OUT "end = $newend\n";
                                                                $sequence = $4;
                                                                $final{$newstart} = $newend;
                                                            }
                                                        }
                                                        my @s = sort {$a <=> $b} keys %final;
                                                        #print OUT "@s\n";
                                                        my @e = sort {$b <=> $a} values %final;
                                                        #print OUT "@e\n";
                                                        $chr=~ s/\-.*$//;
                                                        print OUT qq!$chr\t$s[0]\t$e[0]\t$filename\t1000\t+\t$s[0]\t$e[0]\t0\t$exons{$group}\t$length{$group}\t$offset{$group}\n!;
                                                }	
						else
						{
							if ($strand =~/Plus\s\/\sMinus/)
							{
								#print "$strand\n";
								my $g = 1;
								my %length = ();
								my %offset = ();
								my %exons = ();
								my %store= ();
								my %seq_length = ();
								my $i=0;
								my $r = $i+1;
								my $stop = 1;
								#print "segment count\t";
								#print $segment_count{$chr_};
								while ($i<=$segment_count{$chr_})
								{
									while ($stop == 1)
									{
										$store{qq!$segment_start_s{"$chr_#$i"}#$segment_end_s{"$chr_#$i"}#$segment_seq_s{"$chr_#$i"}!}= $g;
										my $l=$results{"$chr_#$i#end"}-$results{"$chr_#$i#start"}+1;
										my $o=$results{"$chr_#$i#start"}-$results{"$chr_#0#start"};
										my $len = length $segment_seq_q{"$chr_#$i"};
										$seq_length{$g} += $len;
										$length{$g}.="$l,";
										$offset{$g}.="$o,";
										$exons{$g} += 1;
										if ($segment_end_s{"$chr_#$i"} < $segment_start_s{"$chr_#$r"} and $segment_end_q{"$chr_#$i"} > $segment_start_q{"$chr_#$r"})
										{
											$i++;
											$r++;
										}else{
											$stop = 0;
										}
									}
									
									$i++;
									$r++;
									$g++;
									$stop = 1;
								}
								#print OUT "$i, $r, $g\n";
								my @store=();
								foreach my $key (sort {$a <=> $b} keys %store)
								{
								    push @store, "$store{$key}\.$key\n";
								}
								#print OUT "Store:\n@store";
								
								my @length =();
								foreach my $key (sort {$seq_length{$b} <=> $seq_length{$a}} keys %length)
								{
								    push @length, "$key\t$seq_length{$key}\n";
								}
								#print OUT "Length:\n@length";
								my $line = $length[0];
								my $group ="";
								my $length="";
								if ($line =~/^(\d+)\t(\d+)/)
								{
								    $group = $1;
								    #print OUT "Group:$group\n";
								    $length = $2;
								}
								my $newstart = "";
								my $newend = "";
								my $sequence = "";
								my %final = ();
								my $p ="";
								for ($p = 0; $p<$i; $p++)
								{
								    #print OUT "$p\n";
								    if ($store[$p]=~/^($group)\.(\d+)#(\d+)#(.+)$/)
								    {
									$newstart = $2;
									#print OUT "start = $newstart\n";
									$newend = $3;
									#print OUT "end = $newend\n";
									$sequence = $4;
									$final{$newstart} = $newend;
								    }
								}
								my @s = sort {$a <=> $b} keys %final;
								#print OUT "@s\n";
								my @e = sort {$b <=> $a} values %final;
								#print OUT "@e\n";
								$chr=~ s/\-.*$//;
								print OUT qq!$chr\t$s[0]\t$e[0]\t$filename\t1000\t+\t$s[0]\t$e[0]\t0\t$exons{$group}\t$length{$group}\t$offset{$group}\n!;
							}	
						}
					}
				}
				close(OUT);
											
			}
		}
	}
}
