#!/usr/bin/perl
use strict;
##############################################################
# 
# checking variation postition whather change amino acid sequence
# 
# usage : perl check_aa_sequence.pl $ARGV[0] $ARGV[1]  $ARGV[2] 
# 
# perl check_aa_sequence.pl ATCC824_reorganized.fna codon_transfer.txt mutation_in_coding_region.txt 
# perl check_aa_sequence.pl HCCB10218_reorganized.fna codon_transfer.txt mutation_in_coding_region.txt 
# 
# 
# 	$ARGV[0] >>> reference cds position file (***_cds_from_genomic.fna)
# 	$ARGV[1] >>> codon file (codon_transfer.txt)
# 	$ARGV[2] >>> input vcf file (mutation_in_coding_region.txt )
# 
# 
# edit time :
# 20190913
# 20190914
# 20190915
# 20190916
# 
# 

my $header;
my @tmp;
my $tmp;

## mutation file content
my ($chomosome,$variation_position,$type,$Ref,$Alt,$cds_locus_tag,$cds_start_end,$translation_or_not);

my %chomosome;
my %data;

my $tmp_key;

my %position;
my (@pos_tmp, $cds_start_pos, $cds_end_pos);

## mutation modifying related
my ($pos1, $pos2);
my (%cds_way, %cds_mutation_pos, %mutant_edit);
my ($edit_pos, $edit_len);

my @fas_tmp;
my $fasta;
my $fragment;
my %fasta;
my %cds_fasta;

my $locus;

my %tanslation;

my ($aa,$aa_s,$codon,$codon_s,%codon,%codon_abb);

my @cds_mut_pos;
my @mutant_edit;
my $pos_tmp;
my $i;

my ($fragment1, $fragment2, $fragment3);
my $ref_of_alt;

my ($mut_cds_fasta, $cds_fasta, $aa_seq, $mut_aa_seq, $cds_modified_fasta);
my @seq_tmp;

my $file_name;

################################
# 
#  Reading mutation file
# 
open (IN, "$ARGV[2]")||die "$!";

$header=<IN>;
while(<IN>)
{
	chomp;
	@tmp = split "\t",$_;
	
	$chomosome = $tmp[0];
	$variation_position = $tmp[1];
	$type = $tmp[2];
	$Ref = $tmp[3];
	$Alt = $tmp[4];
	$cds_locus_tag = $tmp[5];
	$cds_start_end = $tmp[6];
	$translation_or_not = $tmp[7];
		
	$tmp_key="$chomosome\t$variation_position";
	$data{$tmp_key} = $_;
		
	$locus="$cds_locus_tag\t$cds_start_end";
	


	######################
	# 
	#  Get cds position
	# 
	@pos_tmp=split "-",$cds_start_end;
	
	$cds_start_pos = $pos_tmp[0];
	$cds_end_pos = $pos_tmp[1];
	$pos1 = $cds_start_pos-1;
	
	$tanslation{$locus} = $translation_or_not;
	
	$pos2 = ($cds_end_pos-$cds_start_pos+1);
	
	if ($pos2<0)
	{
		$pos2 = ($cds_start_pos-$cds_end_pos+1);
	}
	
	$cds_fasta{$locus} = "$chomosome-$pos1-$pos2";
	$cds_mutation_pos{$locus} = "$cds_mutation_pos{$locus}---$tmp_key";
	$position{$tmp_key} = $locus;
	
	$edit_pos = ($variation_position-$pos1-1);
	$edit_len = length($Ref);
	$mutant_edit{$tmp_key} = "$cds_locus_tag\t$edit_pos\t$edit_len\t$Ref\t$Alt\t$type";
}

close IN;



######################
# 
# Codon convertion
# 
open (IN,"$ARGV[1]")||die "$!";

while(<IN>)
{
	chomp;
	@tmp = split "\t",$_;
	$codon = $tmp[0];
	$codon_s = lc $codon;
	#print "$codon\t$codon_s\n";
	
	$aa = $tmp[1];
	$aa_s = $tmp[2];
	
	$codon{$codon} = $aa;
	$codon{$codon_s} = $aa;
	
	$codon_abb{$codon} = $aa_s;
	$codon_abb{$codon_s} = $aa_s;
	
	
}

close IN;




#########################
# 
#  Get cds position
# 

system '
	
	if [ ! -d "Sequence_comparison" ] ;then
		mkdir Sequence_comparison
	fi
	cd Sequence_comparison
	
	if [ ! -d "DNA" ] ;then
		mkdir DNA
	fi
	
	if [ ! -d "protein" ] ;then
		mkdir protein
	fi
';


open (CDS, ">./Sequence_comparison/mutation_cds.fna")||die "$!";
open (MutCDS, ">./Sequence_comparison/mutation_cds_modified.fna")||die "$!";

open (faa, ">./Sequence_comparison/mutation_cds.faa")||die "$!";
open (Mutfaa, ">./Sequence_comparison/mutation_cds_modified.faa")||die "$!";

open (Check, ">./Sequence_comparison/amino_acid_primary_check.txt")||die "$!";






foreach $locus (sort keys %cds_fasta)
{
	#################################
	# 
	#  Extraction cds DNA sequence
	# 
	@tmp=split "-",$cds_fasta{$locus};
	$chomosome = $tmp[0];
	$pos1 = $tmp[1];
	$pos2 = $tmp[2];
	
	open (fna, "$ARGV[0]")||die "$!";
	while (<fna>)
	{
		chomp;
		@fas_tmp=split "\t", $_;
		if ($fas_tmp[0] =~ /$chomosome/)
		{
			$fasta = substr $fas_tmp[1], $pos1, $pos2;
		}
	}

	$cds_fasta = $fasta;
	if ( $fasta =~/^ATG/ || $fasta =~/^GTG/ || $fasta =~/^TTG/  || $fasta =~/^ATT/  || $fasta =~/^CTG/  )
	{}  
	elsif ( $fasta =~/^TTA/ || $fasta=~/^CTA/ || $fasta =~/^TCA/ )
	{
		$cds_fasta = reverse $fasta;
		$cds_fasta =~ tr/ATGCatgc/TACGtacg/;
	}
	print CDS ">$locus\n$cds_fasta\n";
	close fna;
	
	
	if ( $tanslation{$locus} =~ /no/ ){}
	elsif ( $tanslation{$locus} =~ /yes/ )
	{
		
		$aa_seq = $cds_fasta;
		$aa_seq =~ s/(...)/"$codon_abb{$1}" || "?"/eg;
	}
	print faa ">$locus\n$aa_seq\n";
	

	
	
	
	#################################
	# 
	#  Mutation cds DNA sequence
	# 
	@cds_mut_pos=split "---", $cds_mutation_pos{$locus}; ## [chromosome]\t[position]
	$fragment = $fasta;
	$tmp = 0;
	for ( $i=1 ; $i < @cds_mut_pos ; $i++ )
	{
		$tmp_key = $cds_mut_pos[$i];
		@mutant_edit = split "\t",$mutant_edit{$tmp_key};
		$pos1 = $mutant_edit[1]-$tmp;		## cutted from postion
		$pos2 = $mutant_edit[2];   			## cutted length
		$Ref = $mutant_edit[3];
		$Alt = $mutant_edit[4];
		$type = $mutant_edit[5];
		
		if ($type=~/SNP/)
		{
			$pos2 = $mutant_edit[2];

			
		}
		elsif(  $type=~/insertion/)
		{
			if(length($Ref)<2)
			{
				$pos2=$mutant_edit[2]-1;
			}
			else
			{
				$tmp = $tmp + ( length($Ref)-length($Alt) );
			}
		}
		elsif( $type=~/deletion/ )
		{
			if(length($Ref)<2)
			{
				$pos2 = $mutant_edit[2]-1;
			}
			else
			{
				$tmp = $tmp - ( length($Ref)-length($Alt) );
			}
		}
		
		$ref_of_alt= substr $fasta, $pos1, $pos2;
		
		if ($ref_of_alt == $Ref )
		{
			$fragment1 = substr $fasta, 0, ($pos1);
			$fragment2 = $Alt;
			$fragment3 = substr $fasta, ($pos1+$pos2), (length($fasta)-$pos1);
			
			$fasta=join "", $fragment1, $fragment2, $fragment3;
			if ($tmp_key =~/4971308/)
			{
				
				#print "gene length: ".length($fasta)."\n";
				#print "position: $pos1\n";
				#print "length:   $pos2\n";
				#print "fragment1 : fasta, 0, $pos1\n";
				#print "fragment2 : fasta, $pos1, $pos2 >> $ref_of_alt\n";
				#print "fragment3 : fasta, ".($pos1+$pos2).", ".(length($fasta)-$pos1)."\n";
			}
			#print MutCDS "$tmp_key\t";
			#print MutCDS "$fasta\n";
		}

	}
	
	$mut_cds_fasta = $fasta;
	if ( $fasta=~/^ATG/ || $fasta=~/^GTG/ || $fasta=~/^TTG/  || $fasta=~/^ATT/  || $fasta=~/^CTG/  )
	{}  
	elsif ( $fasta=~/^TTA/ || $fasta=~/^CTA/ || $fasta=~/^TCA/ )
	{
		$mut_cds_fasta = reverse $fasta;
		$mut_cds_fasta =~ tr/ATGCatgc/TACGtacg/;
	}
	print MutCDS ">$locus\n$mut_cds_fasta\n";
	
	if ( $tanslation{$locus} =~ /no/ ){}
	elsif ( $tanslation{$locus} =~ /yes/ )
	{
		$mut_aa_seq=$mut_cds_fasta;
		$mut_aa_seq =~ s/(...)/"$codon_abb{$1}" || "?"/eg;
	}
	print Mutfaa ">$locus\n$mut_aa_seq\n";
	
	
	
	if ( $mut_aa_seq eq $aa_seq )
	{
		print Check "$locus\tnormal\n";
	}
	else
	{
		print Check "$locus\tamino acid change\n";
		
		$file_name=$locus;
		$file_name=~s/\t/_/g;
		
		open (DNA,">./Sequence_comparison/DNA/$file_name.fa")||die "$!";
		open (protein,">./Sequence_comparison/protein/$file_name.aa")||die "$!";
		
		print DNA ">$locus\_ref\n$cds_fasta\n";
		print DNA ">$locus\_alt\n$mut_cds_fasta\n";
		
		print protein ">$locus\_ref\n$aa_seq\n";
		print protein ">$locus\_alt\n$mut_aa_seq\n";
		
		close DNA;
		close protein;
	}
	
	
	
}

close CDS;
close faa;
close MutCDS;
close Mutfaa;
close Check;









