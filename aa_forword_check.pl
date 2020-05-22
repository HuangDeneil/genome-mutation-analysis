#!/usr/bin/perl
use strict;
##############################################################
# 
# checking variation postition whather change amino acid sequence
# 
# usage : perl aa_forword_check.pl $ARGV[0] $ARGV[1]  $ARGV[2] 
# 
# perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt
# 
# 
# 	$ARGV[0] >>> amino_acid_primary_check.txt
# 	$ARGV[1] >>> mutation reference amino acid sequence (mutation_cds.faa)
# 	$ARGV[2] >>> mutated amino acid sequence (mutation_cds_modified.faa)
# 	$ARGV[2] >>> mutation cds info (mutation_in_coding_region.txt)
# 
# 
# edit time :
# 20190916
# 
# 



my @tmp;
my $tmp;
my $i;

my %list;
my %pseudogenes;
my %ref;
my %alt;
my $count=0;
my $locus;
my @seq;
my $seq;
my $seq_count;
my %transcript_dead;
my $header;
my %aa_change_info;
my $test;


open (IN,"$ARGV[0]") || die "$!";

while(<IN>)
{
	chomp;
	@tmp = split "\t",$_;
	$locus="$tmp[0]\t$tmp[1]";
	$list{$locus}=$tmp[2];
	
	if ($tmp[2]=~/normal/)
	{
		$aa_change_info{$locus}="normal";
		
	}
	elsif ($tmp[2]=~/amino acid change/)
	{
		$aa_change_info{$locus}="amino acid change";
	}
	
}
close IN;




#######################################
# 
#  Reading ref faa
# 
open (IN,"$ARGV[1]") || die "$!";

while(<IN>)
{
	chomp;
	if ($_=~/^>(.+)/)
	{
		$locus=$1;
		$count = 0;
	}
	
	if ($count == 0)
	{
		
		$count++;
	}
	else
	{
		$seq=$_;
		$ref{$locus}=$seq;
		$seq_count=0;
		@seq=split "",$seq;
		
		for ($i=0;$i<@seq;$i++)
		{
			if ( $seq[$i] =~ /\*/ )
			{
				$seq_count++;
			}
		}
		
		if ( $seq_count > 1 || $seq_count == 0 )
		{
			$pseudogenes{$locus}="origin_sequence_dead";
			
		}
		else
		{
			$pseudogenes{$locus}="no";
		}
	}
}
close IN;

# foreach $locus(sort keys %ref)
# {
	# print "$locus\t$pseudogenes{$locus}\n$ref{$locus}\n";
	
# }

#######################################
# 
#  Reading alt faa
# 
open (IN,"$ARGV[2]") || die "$!";

while(<IN>)
{
	chomp;
	if ($_=~/^>(.+)/)
	{
		$locus=$1;
		$count = 0;
		#print"$locus\n";
	}
	
	if ($count == 0)
	{
		$count++;
	}
	else
	{
		$seq=$_;
		$alt{$locus}=$seq;
		$seq_count=0;
		@seq=split "",$seq;
		
		for ($i=0;$i<@seq;$i++)
		{
			if ( $seq[$i] =~ /\*/ )
			{
				$seq_count++;
			}
		}
		
		
		if ( $pseudogenes{$locus} eq "origin_sequence_dead" )
		{
			
			if ($seq_count == 1)
			{
				if($seq=~/[A-Z]+\*$/){
					
					$pseudogenes{$locus}= "sequence_relived";
				}
			}
			$transcript_dead{$locus}="no";
		}
		elsif ( $seq_count > 1 || $seq_count == 0 )
		{
			$transcript_dead{$locus} = "dead";
		}
		else
		{$transcript_dead{$locus}="no";}
	}
}
close IN;






foreach $locus (sort keys %list)
{
	if ($list{$locus} eq "amino acid change")
	{
		if ($aa_change_info{$locus} eq "normal" )
		{}
		elsif ($pseudogenes{$locus} eq "origin_sequence_dead" )
		{
			$aa_change_info{$locus}="Maybe is a pseudogene";
		}
		elsif ($transcript_dead{$locus} eq "dead")
		{
			$aa_change_info{$locus}="multiple or lose stop codon";
		}
		elsif ($pseudogenes{$locus} eq "sequence_relived" )
		{
			$aa_change_info{$locus}="aa changed but becomed a translatable sequece";
		}
		else
		{
			$aa_change_info{$locus}="amino acid change but maybe work";
		}
	}
}
#print "$test\n";




#################################
# 
#  Reading mutation cds info
# 

my %data;
my %mut_count;
my $report;

open (IN,"$ARGV[3]") || die "$!";
$header=<IN>;
while(<IN>)
{
	chomp;
	@tmp=split "\t",$_;
	$locus="$tmp[5]\t$tmp[6]";
	
	if ($pseudogenes{$locus} eq "origin_sequence_dead" )
	{
		if ($tmp[9]=~/pseudogene/)
		{
			$aa_change_info{$locus}="Oriental recognized as pseudogene and dead gene";
		}
	}
	$mut_count{$locus}++;
	$report=$aa_change_info{$locus};
	$report=~s/\t/ in /;
	$data{$locus}="$tmp[5]\t$tmp[6]\t$mut_count{$locus}\t$tmp[7]\t$report\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\t$tmp[13]\t$tmp[14]";
	
}
close IN;

open (OUT,">final_report.txt")||die"$!";
$header="cds_locus_tag\tcds_start-end\tmutation_point\ttranslation\taa_change_info\tID\tproduct name\tproduct name(uniprot)\tEC\tKEGG (KO)\tGene ontology (biological process)\tGene ontology (cellular component)\tGene ontology (molecular function)";

print OUT "$header\n";
foreach (sort keys %data)
{
	print OUT "$data{$_}\n"
}
close OUT;




