#!/usr/bin/perl -w
use strict;

my %pep = ();
my $pep="";
my $acc="";
my $dir="";
my $file="";
my $mouse_uni="";
my $rat_uni="";
my $mouse_prot="";
my $mouse_prot_="";
my $rat_prot="";
my $rat_prot_="";
my $mouse_tran="";
my $rat_tran="";
my %mus_uni_to_prot=();
my %rat_uni_to_prot=();
my %mus_prot_to_tran=();
my %rat_prot_to_tran=();
my @temp = ();
my @ensm_p=();
my @ensm_t=();
my $uni = "";
my %uni = ();
my %ensm_mouse = ();
my %ensm_rat = ();
my %ensm_p=();
my $tran = "";
my $seq = "";
my %seq = ();
my @seq = ();

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }
if ($ARGV[1]=~/\w/) { $file=$ARGV[1];} else { print "no file"; }
if (opendir(DIR,"$dir")){
        my @allfiles=readdir DIR;
        closedir DIR;
        foreach my $filename (@allfiles){
		if ($filename=~/\.out$/i){
			print qq!$filename\n!;
			open (IN_,"$dir/$filename") || die "Could not open $dir/$filename\n";
			while (my $line = <IN_>){
			    if ($line =~ /^(.+)\s+(.+)$/){
				$pep = $1;
				#print "$pep\n";
				$acc = $2;
                                #print "$acc\n";
				$pep{$pep}=$acc;
			    }
			}
		}
	}
}
open (IN, "$dir/$file");
while (my $line = <IN>){
	chomp $line;
	#print "$line\n";
	foreach (keys %pep){
	#print "peptides: $_\n";	
		if ($line eq $_) {
			#print "$line equals $_\n";
			push @temp, $pep{$line};
		}
	}
}
#foreach my $line (@temp){
#    chomp $line;
#    $line =~ s/#/\n/g;

#}

#print "temp:\n @temp";
foreach my $line (@temp){
	#print "temp line:$line\n";
	chomp $line;
	if ($line=~/^(.{6})$/){
		$uni = $1;
		#print "uni: $uni\n";
		$uni{$uni}=1;
		#print "from temp:$uni\n";
	}elsif ($line=~/^(ENSMUSP.+)$/){
		my $em = $1;
		#print "$em\n";
		$ensm_mouse{$em} = 1;
		#print "ensm mouse from temp:$em\n";
	}#elsif ($line=~/^(ENSRNOP.+)/){
		#my $er = $1;
		#$ensm_rat{$er}=1;
		#print "ensm ratfrom temp:$\n";
	#}
}

if (opendir(DIR,"$dir")){
        my @allfiles=readdir DIR;
        closedir DIR;
        foreach my $filename (@allfiles){
		if ($filename=~/\.txt$/i){
		    #print qq!$filename\n!;
		    open (IN_,"$dir/$filename") || die "Could not open $dir/$filename\n";
		    while (my $line = <IN_>){
			    if ($line=~/^(.{6})\s+(ENSMUSP.+)$/){
				    $mouse_uni = $1;
				    #print "mouse uni; $mouse_uni\n";
				    $mouse_prot = $2;
				    $mus_uni_to_prot{$mouse_uni}=$mouse_prot;
				    #print "uni to ensm mouse works:$mouse_uni:$mouse_prot\n";
				    if (exists $uni{$mouse_uni}){
					    push @ensm_p, $mus_uni_to_prot{$mouse_uni};
					    #print "found mouse match\n";
					}
			    }elsif ($line=~/^(.{6})\s+(ENSRNOP.+)$/){
				    $rat_uni = $1;
				    #print "rat uni:$rat_uni\n";
				    $rat_prot = $2;
				    $rat_uni_to_prot{$rat_uni}=$rat_prot;
				    #print "uni to ensm rat works:$rat_uni:$rat_prot\n";
				    if(exists $uni{$rat_uni}){
					    push @ensm_p, $rat_uni_to_prot{$rat_uni};
					    #print "found rat match\n";
					}
			    }elsif ($line=~/^(ENSRNOP.+)\s+(ENSRNOT.+)$/){
				    $rat_prot_ = $1;
				    #print "rat ensm prot:$rat_prot_\n";
				    $rat_tran = $2;
				    $rat_prot_to_tran{$rat_prot_}=$rat_tran;
				    #print "prot and tran storage in rat works:$rat_prot_:$rat_tran\n";
			    }elsif ($line =~/^(ENSMUSP.+)\s+(ENSMUST.+)$/){
				    $mouse_prot_ = $1;
				    #print "mouse ensm prot: $mouse_prot_\n";
				    $mouse_tran = $2;
				    $mus_prot_to_tran{$mouse_prot_}=$mouse_tran;
				    #print "prot and tran storage in mouse works:$mouse_prot_:$mouse_tran\n";
				}
			}
		}
	}
} else {print "could not open $dir\n"};


foreach my $key (@ensm_p){
	#print "$key\n";
	$ensm_p{$key}=1;
}
foreach (keys %ensm_p){
	chomp $_;
	#print "key: $_\n";
	if (exists $rat_prot_to_tran{$_}){
		push @ensm_t, $rat_prot_to_tran{$_};
		#print "ensm prot to trans from uni in rat works:$_\n";
	}elsif (exists $mus_prot_to_tran{$_}){
		push @ensm_t, $mus_prot_to_tran{$_};
		#print "ensm prot to trans from uni in mouse works:$_\n";
	}
}
foreach (keys %ensm_mouse){
	if (exists $mus_prot_to_tran{$_}){
		push @ensm_t, $mus_prot_to_tran{$_};
		#print "ensm prot to trans in mouse works:$_\n";
	}
}
#foreach (keys %ensm_rat){
	#if (exists $rat_prot_to_tran{$_}){
		#push @ensm_t, $rat_prot_to_tran{$_};
		#print "ensm prot to trans in rat works:$_\n";
	#}
#}

open (RSEQ, "$dir/rat_ensm_with_seq.txt");
open (MSEQ, "$dir/mouse_ensm_with_seq.txt");

while (my $line = <RSEQ>){
	$line =~/^\s*>(.+)\s+(\w+)$/;
	$tran = $1;
	#print "rat transcript $tran\n";
	$seq = $2;
	$seq{$tran}=$seq;
	#print "seq storage in rat works:$seq{$tran}\n";
}
while (my $line = <MSEQ>){
	$line =~/^\s*>(.+)\s+(\w+)$/;
	$tran = $1;
	#print "mouse transcript $tran\n";
	$seq = $2;
	$seq{$tran}=$seq;
	#print "seq storage in mouse works:$seq{$tran}\n";
}

foreach my $id (@ensm_t){
	#print "$id\n";
	if (exists $seq{$id}){
		push @seq, ">$id\t$seq{$id}\n";
		#print "final doc creation works\n";
	}
}
open (OUT, ">>$dir/sequences.txt");
print OUT @seq;
