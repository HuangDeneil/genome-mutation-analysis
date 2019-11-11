#!/usr/bin/perl
use strict;
##############################################################
# 
# This perl is summary of mutation position infrot of genes
# 
# usage : perl noncoding_forward_analysis.pl $ARGV[0] 
# 
# perl noncoding_forward_analysis.pl mutation_in_non-coding_region_get_nearby.txt
# 
# 
# 	$ARGV[0] >>> noncoding_analysis.pl result file (mutation_in_non-coding_region_get_nearby.txt)
# 
# 
# mutation_in_non-coding_region_get_nearby.txt
# ----------------------------------
# Chomosome	variation position	variation type	REF	ALT	location	nearby_cds1	nearby_cds2	nearby_cds1_info	nearby_cds2_info
# NC_003030.1	1089509	insertion	TAA	TAAA	non-coding region region	CA_C0950-head	CA_C0949-tail	hypothetical protein	transcriptional regulator
# NC_003030.1	1110036	deletion	AAG	A	non-coding region region	CA_Cr029-tail	CA_Cr030-head	23Sj ribosomal RNA	5Sj ribosomal RNA
# NC_003030.1	1110047	SNP	A	G	non-coding region region	CA_Cr029-tail	CA_Cr030-head	23Sj ribosomal RNA	5Sj ribosomal RNA
# NC_003030.1	1228591	SNP	G	T	non-coding region region	CA_C1075-tail	CA_C1076-head	beta-glucosidase	hypothetical protein
# NC_003030.1	1395215	deletion	GAAA	GAA	non-coding region region	CA_C1249-tail	CA_C1251-head	septum site-determining protein MinD%2C ATPase	cell cycle protein FtsW
# NC_003030.1	1395357	SNP	G	A	non-coding region region	CA_C1249-tail	CA_C1251-head	septum site-determining protein MinD%2C ATPase	cell cycle protein FtsW
# 
# -----------------------------------
# 
# 
# edit time :
# 20190917
# 
# 

my $header;
my $key_tmp;
my @tmp;
my %data;
my %cds_info;
my %count;
my $locus;
my $loc1;
my $loc2;

open(IN,"$ARGV[0]")||die "$!";

while (<IN>){
	chomp;
	@tmp=split "\t",$_;
	
	$key_tmp="$tmp[0]-$tmp[1]";
	
	$data{$key_tmp}=$_;
	
	if ($tmp[6]=~/head/){
		
		$locus="$tmp[6]";
		$locus=~s/-head//;
		$count{$locus}++;
		
		$cds_info{$locus}=$tmp[8];
		
		if($tmp[7]=~/head/)
		{
			$locus="$tmp[7]";
			$locus=~s/-head//;
			$count{$locus}++;
			$cds_info{$locus}=$tmp[9];
		}
	}
	elsif($tmp[7]=~/head/)
	{
		$locus="$tmp[7]";
		$locus=~s/-head//;
		$count{$locus}++;
		$cds_info{$locus}=$tmp[9];
	}
}

close IN;

print "locus tag\tmutation in front of locus\tlocus product\n";

foreach $locus (sort keys %count)
{
	print "$locus\t$count{$locus}\t$cds_info{$locus}\n";
}







