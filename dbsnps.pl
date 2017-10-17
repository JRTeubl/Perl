#!/usr/bin/perl -w
use strict;

my $dir = "";
my $whim_germ="";
my $whim_som="";
my $dbsnp = "";
my $whim="";
my $germline = "germline";
my $somatic = "somatic";
my $chr_whim_germ = "";
my $chr_whim_som = "";
my $chr_whim ="";
my $pos_whim_germ = "";
my $pos_whim_som = "";
my $chr_dbsnp = "";
my $pos_dbsnp = "";
my $chr = "";
my %hash_whim = ();
my %hash_dbsnp= ();
my $whim_pos = "";
my $whim_dna = "";
my $pos = "";
my $line = "";

if ($ARGV[0]=~/\w/) {$dir=$ARGV[0]};
if ($ARGV[1]=~/\w/) {$dbsnp=$ARGV[1]};
if ($ARGV[2]=~/\w/) {$whim_germ=$ARGV[2]};
if ($ARGV[3]=~/\w/) {$whim_som=$ARGV[3]};

if ($whim_germ=~/^(whim\d+)\.(germline)\.snvs.vcf$/){
    $whim = $1;
    #print "$whim_germ\n";
}

open (IN, $whim_germ);

while (my $line_ = <IN>){
    if ($line_=~ /^(\d+)\s+(\d+).*$/){
        $chr_whim_germ = $1;
        $pos_whim_germ = $2;
        $pos_whim_germ.= $germline;
        #print "$pos_whim_germ\n";
        $hash_whim{$pos_whim_germ}=$chr_whim_germ;
    } 
}

open (IN1, $whim_som);
while (my $line_ = <IN1>){
    if ($line_=~ /^(\d+)\t(\d+).*$/){
        $chr_whim_som = $1;
        $pos_whim_som = $2;
        $pos_whim_som.= $somatic;
        #print "$pos_whim_som\n";
        $hash_whim{$pos_whim_som}=$chr_whim_som;
    }
}

open (IN2, $dbsnp);
while (my $line_ = <IN2>){
    if ($line_=~/^(\d+)\t(chr\d+)\t(\d+)\t(.*\t){20}(.*)$/){
        $chr_dbsnp = $2;
        $pos_dbsnp = $3;
        #print "$pos_dbsnp\n";
        if ($chr_dbsnp =~ /^(chr)(\d+)$/){
            $chr = $2;
            $hash_dbsnp{$pos_dbsnp}=$chr_dbsnp
        }
    }
}

foreach my $keydb (keys %hash_dbsnp){
    print "db:$keydb\n";
    foreach my $keywhim (keys %hash_whim){
        if ($keywhim =~ /^(\d+)(\w+)$/){
            $whim_pos = $1;
            print "whim:$whim_pos\n";
            $whim_dna = $2;
        }
        if ($keydb = $whim_pos){
            delete $hash_whim{$keydb};
        }
    }
}

open (OUT, ">>$dir/$whim.new.txt");

foreach my $key (keys %hash_whim){
    if ($key =~/^(\d+)(\w+)$/){
        $pos = $1;
        $line = $2;
        if ($line =~/germline/){
            print OUT "$hash_whim{$key}\t$pos\t0\n";
        }elsif ($line =~ /somatic/){
            print OUT "$hash_whim{$key}\t$pos\t1\n";
        }
    }
}