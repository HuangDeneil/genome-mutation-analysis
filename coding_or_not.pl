#!/usr/bin/perl
use strict;
#############################################
# 
# This perl script split coding and non-coding mutation
# 
# usage : perl coding_or_not.pl $ARGV[0] 
# 
# $ARGV[0] >>> mutation cds file
# 
# perl coding_or_not.pl cds_merged.txt
# 
# 
# Edit time 
# 20190910
# 
# 
# 


my $tmp;
my @tmp;
my $cds_name;
my %coding;
my %noncoding;
my $header;

##############################
# 
# Doing uniprot annotation
# 
open (IN,"$ARGV[0]")||die "$!";
$header=<IN>;
while(<IN>)
{
	chomp;
	@tmp= split "\t",$_;
	
	if ($tmp[1]==""){}
	elsif($tmp[0]=~/^#/){}
	elsif ($tmp[5]=~/non-coding region/)
	{
		
		$noncoding{"$tmp[0]\t$tmp[1]"}=$_;
		
	}
	else
	{
		$coding{"$tmp[0]\t$tmp[1]"}=$_;
	}
	

}


close IN;

open (OUTcoding,">mutation_in_coding_region.txt")||die"$!";
print OUTcoding "$header";
foreach (sort keys %coding)
{
	print OUTcoding "$coding{$_}\n";
	
}

close OUTcoding;

open (OUTnoncoding,">mutation_in_non-coding_region.txt")||die"$!";
print OUTnoncoding "$header";
foreach (sort keys %noncoding)
{
	print OUTnoncoding "$noncoding{$_}\n";
	
}

close OUTnoncoding;
