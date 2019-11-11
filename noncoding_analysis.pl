#!/usr/bin/perl
use strict;
##############################################################
# 
# This perl is finding nearby cds of mutation position
# 
# usage : perl noncoding_analysis.pl $ARGV[0] $ARGV[1] 
# 
# perl noncoding_analysis.pl GCF_000008765.1_ASM876v1_genomic.gff mutation_in_non-coding_region.txt
# 
# 
# 	$ARGV[0] >>> reference cds position file (***_cds_from_genomic.fna)
# 	$ARGV[1] >>> input noncoding mutation file (mutation_in_non-coding_region.txt)
# 
# 
# mutation_in_non-coding_region.txt:
# ----------------------------------
# Chomosome	variation position	variation type	REF	ALT cds_locus_tag
# NC_003030.1	1089509	insertion	TAA	TAAA	non-coding region
# NC_003030.1	1110036	deletion	AAG	A	non-coding region
# NC_003030.1	1110047	SNP	A	G	non-coding region
# 
# -----------------------------------
# 
# 
# edit time :
# 20190913
# 
# 

#####################################################################
# 
# Print current time 
# 
# system'
	# now=`date "+%Y-%m-%d_%H:%M:%S"`
	# printf "%s\n" "$now"
	# ';

if (!$ARGV[0]||!$ARGV[1]){
	print "\nUsage error!!\n\n";
	print "usage : perl noncoding_analysis.pl \$ARGV[0] \$ARGV[1] \n\n";
	print "\$ARGV[0] >>> reference cds position file (GCF_000008765.1_ASM876v1_genomic.gff)\n";
	print "\$ARGV[1] >>> nput noncoding mutation file (mutation_in_non-coding_region.txt)\n";
	die"\n";
	
}



my @tmp;
my ($ref_len,$alt_len);
my $chrmosome;
my %chrmosome;				## use chromosome_position as key stroge chromosome ID
my %noncoding_input_data; 	## use chromosome_position as key stroge mutation info
my %locus_tag;		# record cds position
my $locus;
my $position;
my $position1;
my $position2;
my $key_tmp;
my $header;
my $info;



###############################
# 
# read noncoding mutation position
# 
open (INnoncoding,"$ARGV[1]")||die "$!";

chomp ($header=<INnoncoding>);
while (<INnoncoding>)
{
	chomp;
	
	if (/^\#/){}
	else
	{
		@tmp= split "\t",$_;
		$key_tmp="$tmp[0]_$tmp[1]";				### chromosome_mutation-position
		$chrmosome{$key_tmp}=$tmp[2];
		$info="$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]";
		$noncoding_input_data{$key_tmp}=$info;
		
	}
	
}

close INnoncoding;





###############################
# 
# genomic cds position
# 
open (INcds,"$ARGV[0]")||die "$!";
my $nearby1;
my $nearby2;
my %nearby1;
my %nearby2;
my $nearby_tmp;
my @nearby_tmp;
my %nearby;

my $values_start;
my $values_end;
my $values_last1;
my $values_last2;

my @tmp_mut;
my $cds_start;
my $cds_end;
my %cds_pos_to_locus;
my %cds_pos_to_locus_pos;
my %gene_id;
my $gene_id;
my %protein_ID;
my $protein_ID;
my $protein_value;



###################
## get cds position
while(<INcds>)
{
	chomp;
	
	if (/^\#/){}
	else{
		@tmp= split "\t",$_;
		
		$chrmosome=$tmp[0];
		
		if ($tmp[2]=~/region/)
		{
			$locus=$tmp[0];
			$position="$tmp[3]-$tmp[4]";
			$locus_tag{$locus}=$position;
			$cds_start=$tmp[3];
			$cds_end=$tmp[4];
			
			foreach $key_tmp(sort keys %noncoding_input_data)
			{
				if ($key_tmp=~/$tmp[0]/)
				{
					
					@tmp_mut=split "\t",$noncoding_input_data{$key_tmp};
					
					if ( (1/(($tmp_mut[1]-$cds_start)**2)) > (1/(($tmp_mut[1]-$cds_end)**2)) )
					{
						$nearby1=$cds_start;
						$nearby2=$cds_end;
						$nearby1{$key_tmp}=$nearby1;
						$nearby2{$key_tmp}=$nearby2;
					}
					else
					{
						$nearby1=$cds_end;
						$nearby2=$cds_start;
						$nearby1{$key_tmp}=$nearby1;
						$nearby2{$key_tmp}=$nearby2;
					}
				}
			}
		}
		elsif ( $tmp[2]=~/gene/ )
		{
			
			##############################
			# 
			#  Get cds locus_tag			
			# 
			if ( $tmp[8]=~/\;locus_tag=(.+?)\;/)
			{
				$locus = $1;

			}
			elsif ( $tmp[8]=~/\;locus_tag=(.+)/ )
			{
				$locus = $1;
			}
			
			if ($tmp[8]=~/\;Dbxref=(.+?)\;/)
			{
				$gene_id{$1}=$locus;
			}
			
			if(/\;gene_biotype=pseudogene\;/)
			{
				$protein_value="\tpseudogene";
				$protein_ID{$locus}=$protein_value;
			}
			
			
			
			##############################
			#
			#  Read cds postion
			# 
			$position="$chrmosome-$tmp[3]-$tmp[4]";	## chrmosome-cds_start-cds_end position
			$locus_tag{$locus}=$position;
			$cds_start=$tmp[3];				## cds start position
			$cds_end=$tmp[4];				## cds end position
			
			$cds_pos_to_locus{"$chrmosome-$tmp[3]"}="$locus-head";
			$cds_pos_to_locus{"$chrmosome-$tmp[4]"}="$locus-tail";
			
			
			
			
			#############################
			#
			#  Calculating distance
			#
			foreach $key_tmp(sort keys %noncoding_input_data)
			{
				
				@tmp_mut=split "\t",$noncoding_input_data{$key_tmp};
				
				$values_start=($cds_start-$tmp_mut[1]);
				$values_end=($cds_end-$tmp_mut[1]);
				
				if ( $values_start < 0 || $values_end < 0 )
				{
					$values_start=($tmp_mut[1]-$cds_start);
					$values_end=($tmp_mut[1]-$cds_end);
				}
				
				
				
				
				
				
				#####################################
				#
				#  Find distance under 3000 bp
				#
				if ( $values_start < 3000  || $values_end < 3000 )
				{
					$values_last1=($tmp_mut[1]-$nearby1{$key_tmp});
					$values_last2=($tmp_mut[1]-$nearby2{$key_tmp});
					
					if ( $values_last1 < 0 )
					{
						$values_last1=($nearby1{$key_tmp}-$tmp_mut[1]);
					}
					if ( $values_last2 < 0 )
					{
						$values_last2=($nearby2{$key_tmp}-$tmp_mut[1]);
					}
					
					
					if ( $values_start < $values_end )
					{		### check which cds start or end position more nearby with mutation 
						
						if ($values_start < $values_last1 )
						{	### check which cds start or end position more nearby with
						
							$nearby2=$nearby1{$key_tmp};
							$nearby1=$cds_start;
							$nearby1{$key_tmp}=$nearby1;
							$nearby2{$key_tmp}=$nearby2;
						}
						elsif ($values_start < $values_last2)
						{
							$nearby2=$cds_start;
							$nearby2{$key_tmp}=$nearby2;
						}
					}
					elsif ( $values_start > $values_end )
					{		
						if ($values_end < $values_last1 )
						{	
							$nearby2=$nearby1{$key_tmp};
							$nearby1=$cds_end;
							$nearby1{$key_tmp}=$nearby1;
							$nearby2{$key_tmp}=$nearby2;
						}
						elsif ($values_end < $values_last2)
						{
							$nearby2=$cds_end;
							$nearby2{$key_tmp}=$nearby2;
						}
					}
					$nearby_tmp="$chrmosome-$nearby1{$key_tmp}\t$chrmosome-$nearby2{$key_tmp}";
					$nearby{$key_tmp}=$nearby_tmp;
				}
			}
		}
		elsif($tmp[2]=~/CDS/ || $tmp[2]=~/RNA/ )
		{
			#################
			#
			# Get CDS info
			#
			if ($tmp[8]=~/\,GeneID:(.+?)\;/)
			{
				$gene_id="GeneID:$1";
				$locus=$gene_id{$gene_id};
			}
			
			if ($tmp[8]=~/\;protein_id=(.+?)\;/)
			{
				$protein_ID=$1;
			}
			
			if($tmp[8]=~/\;product=(.+?)\;/)
			{
				$protein_value=join "\t",$protein_ID,$1;
				$protein_ID{$locus}=$protein_value;
			}
			elsif($tmp[8]=~/\;product=(.+)/)
			{
				$protein_value="\t$1";
				$protein_ID{$locus}=$protein_value;
			}
		}
	}
}
close INcds;





###################################################
# 
# Reorganized nearby cds information & output file
# 

my $pos_tmp;
my @pos_tmp;
my $i;
my @nearby_id;
my $nearby_id;



open (OUT,">mutation_in_non-coding_region_get_nearby.txt")||die"$!";

print OUT "Chomosome\tvariation position\tvariation type\tREF\tALT\t";
print OUT "location\tnearby_cds1\tnearby_cds2\tnearby_cds1_info\tnearby_cds2_info\n";

foreach $key_tmp(sort keys %noncoding_input_data)
{
	@tmp_mut=split "\t",$noncoding_input_data{$key_tmp};
	@nearby_tmp=split "\t",$nearby{$key_tmp};
	
	print OUT "$noncoding_input_data{$key_tmp}";
	
	for($i=0;$i<@nearby_tmp;$i++)
	{
		print OUT "\t$cds_pos_to_locus{$nearby_tmp[$i]}";
	}
	
	for($i=0;$i<@nearby_tmp;$i++)
	{
		@tmp =split "-",$cds_pos_to_locus{$nearby_tmp[$i]};
		$locus=$tmp[0];
		
		@tmp =split "\t",$protein_ID{$locus};
		
		print OUT "\t$tmp[1]";
		
	}
	
	
	
	print OUT "\n";
	
	@nearby_id=qw//;
}


close OUT;






