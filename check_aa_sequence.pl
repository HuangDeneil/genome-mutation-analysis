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
my $chomosome;
my %chomosome;
my %data;
my $tmp_key;
my %position;
my $position;
my $cds_start_pos;
my $cds_end_pos;
my @pos_tmp;
my $pos1;
my $pos2;
my %cds_way;
my %cds_mutation_pos;
my %mutant_edit;
my $edit_pos;
my $edit_len;

my @fas_tmp;
my $fasta;
my $fragment;
my %fasta;
my %cds_fasta;

my $Alt;
my $Ref;
my $type;
my $locus;

my %tanslation;

################################
# 
#  Reading mutation file
# 
open (IN, "$ARGV[2]")||die "$!";

$header=<IN>;
while(<IN>)
{
	chomp;
	@tmp=split "\t",$_;
	$tmp_key="$tmp[0]\t$tmp[1]";
	
	$data{$tmp_key}=$_;
	$chomosome=$tmp[0];
	$type=$tmp[2];
	$Ref=$tmp[3];
	$Alt=$tmp[4];
	$locus="$tmp[5]\t$tmp[6]";
	
	######################
	# 
	#  Get cds position
	# 
	$position=$tmp[6];
	@pos_tmp=split "-",$position;
	
	$cds_start_pos=$pos_tmp[0];
	$cds_end_pos=$pos_tmp[1];
	$pos1=$cds_start_pos-1;
	
	$tanslation{$locus}=$tmp[7];
	
	$pos2=($pos_tmp[1]-$pos_tmp[0]+1);
	
	if ($pos2<0)
	{
		$pos2=($pos_tmp[0]-$pos_tmp[1]+1);
	}
	
	$cds_fasta{$locus}="$chomosome-$pos1-$pos2";
	$cds_mutation_pos{$locus}="$cds_mutation_pos{$locus}---$tmp_key";
	$position{$tmp_key}=$locus;
	
	$edit_pos=($tmp[1]-$pos1-1);
	$edit_len=length($Ref);
	$mutant_edit{$tmp_key}="$tmp[5]\t$edit_pos\t$edit_len\t$Ref\t$Alt\t$type";
}

close IN;

######################
# 
# Codon convertion
# 
my $aa;
my $aa_s;
my $codon;
my $codon_s;
my %codon;
my %codon_abb;



open (IN,"$ARGV[1]")||die "$!";

while(<IN>)
{
	chomp;
	@tmp=split "\t",$_;
	$codon = $tmp[0];
	$codon_s= lc $codon;
	#print "$codon\t$codon_s\n";
	
	$aa = $tmp[1];
	$aa_s = $tmp[2];
	
	$codon{$codon}=$aa;
	$codon{$codon_s}=$aa;
	
	$codon_abb{$codon}=$aa_s;
	$codon_abb{$codon_s}=$aa_s;
	
	
}

close IN;




#########################
# 
#  Get cds position
# 

my @cds_mut_pos;
my @mutant_edit;
my $pos_tmp;
my $i;

my $fragment1;
my $fragment2;
my $fragment3;
my $ref_of_alt;

my $mut_cds_fasta;
my $cds_fasta;
my $aa_seq;
my $mut_aa_seq;
my $cds_modified_fasta;
my @seq_tmp;

my $file_name;

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
open (faa, ">./Sequence_comparison/mutation_cds.faa")||die "$!";
open (MutCDS, ">./Sequence_comparison/mutation_cds_modified.fna")||die "$!";
open (Mutfaa, ">./Sequence_comparison/mutation_cds_modified.faa")||die "$!";
open (Check, ">./Sequence_comparison/amino_acid_primary_check.txt")||die "$!";






foreach $locus (sort keys %cds_fasta)
{
	#################################
	# 
	#  Extraction cds DNA sequence
	# 
	@tmp=split "-",$cds_fasta{$locus};
	$chomosome=$tmp[0];
	$pos1=$tmp[1];
	$pos2=$tmp[2];
	
	open (fna, "$ARGV[0]")||die "$!";
	while (<fna>)
	{
		chomp;
		@fas_tmp=split "\t",$_;
		if ($fas_tmp[0]=~/$chomosome/)
		{
			$fasta=substr $fas_tmp[1], $pos1, $pos2;
		}
	}
	
	$cds_fasta = $fasta;
	if ( $fasta=~/^ATG/ || $fasta=~/^GTG/ || $fasta=~/^TTG/  || $fasta=~/^ATT/  || $fasta=~/^CTG/  )
	{}  
	elsif ( $fasta=~/^TTA/ || $fasta=~/^CTA/ || $fasta=~/^TCA/ )
	{
		$cds_fasta = reverse $fasta;
		$cds_fasta =~ tr/ATGCatgc/TACGtacg/;
	}
	print CDS ">$locus\n$cds_fasta\n";
	close fna;
	
	
	if ( $tanslation{$locus} =~ /no/ ){}
	elsif ( $tanslation{$locus} =~ /yes/ )
	{
		
		$aa_seq=$cds_fasta;
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
	for ($i=1;$i<@cds_mut_pos;$i++)
	{
		$tmp_key=$cds_mut_pos[$i];
		@mutant_edit=split "\t",$mutant_edit{$tmp_key};
		$pos1=$mutant_edit[1]-$tmp;
		$pos2=$mutant_edit[2];
		$Ref=$mutant_edit[3];
		$Alt=$mutant_edit[4];
		$type=$mutant_edit[5];
		
		if ($type=~/SNP/)
		{
			$pos2=$mutant_edit[2];
		}
		elsif(length($Ref)<2)
		{
			$pos2=$mutant_edit[2]-1;
		}
		else
		{
			$tmp = $tmp + ( length($Ref)-length($Alt) );
		}
		
		
		$ref_of_alt= substr $fasta, $pos1, $pos2;
		
		if ($ref_of_alt == $Ref)
		{
			$fragment1= substr $fasta, 0, ($pos1);
			$fragment2= $Alt;
			$fragment3= substr $fasta, ($pos1+$pos2), (length($fasta)-$pos1);
			
			$fasta=join "", $fragment1, $fragment2, $fragment3;
			
			if ($tmp_key=~/4972538/)
			{
				print "$cds_mutation_pos{$locus}\n";
				print "mutant_edit:\n@mutant_edit\n";
				print length($fasta)."\n";
				print "position: $pos1\n";
				print "length:   $pos2\n";
				print "fragment1 : fasta, 0, $pos1\n";
				print "fragment2 : fasta, $pos1, $pos2 >> $ref_of_alt\n";
				print "fragment3 : fasta, ".($pos1+$pos2).", ".(length($fasta)-$pos1)."\n";
			}
			
			print MutCDS "$tmp_key\t";
			print MutCDS "$fragment1\n";
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
		print Check "$locus\tanimo acid change\n";
		
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









