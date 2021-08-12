#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $now;
my $PROG = basename $0;
my $exit_code;
my ($input, $input1, $input2, $input3, $input4, $output, $output1, $output2, $output3, $output4);
my ($ResultPath, $outpath);
my ($gff, $vcf_file);

sub display_version {
  print STDERR " 
  ###############= VERSION =################
    2021/08/10
    $PROG v1.0.0 (beta version), 
    hudeneil (hudeneil\@gmail.com)
   ";
  exit 0;
}

sub usage {
    my $exit_code = @_ ? shift : 64;
    print STDERR <<EOF;
This program is checking variant position (annotation step):

Usage: $PROG [options] 
Options:
Input/Output:
  
  -vcf STR        		  Input variant calling file (vcf) 
  -gff STR         		  Reference genome gff
  
  --Result,-o  STR        output file [ default: variantion_all.txt ]
  
  --help,-h               Print this message
  --version,-v            Print version information

example:

perl $PROG \
-gff genome.gff \
-vcf input.vcf \
-o ./result/test_sample 

EOF
    exit $exit_code;
}

GetOptions(
    "h" => \&display_help,
    "v" => \&display_version,
    "vcf=s" => \$vcf_file,
    "gff=s" => \$gff,
    "o=s" => \$ResultPath,
    
    "help" => \&display_help,
    "version" => \&display_version,
    "Result=s" => \$ResultPath,
);
# Check input exists
if ( ! defined $ResultPath ) { $ResultPath = "variantion_all.txt";}
if ( ! defined $vcf_file ) { print STDERR "\n !!!!! please input variant calling file (vcf) !!!!!! \n"; usage(0);}
if ( ! defined $gff ) { print STDERR "\n !!!!! please input reference gff !!!!!! \n"; usage(0);}


## Check input file exists
if ( ! -f "$vcf_file" ) { print STDERR "\n\n!!!! $vcf_file : vcf file is not exists !!!!!\n\n";exit 1;}
if ( ! -f "$gff" ) { print STDERR "\n\n!!!! $gff : gff file is not exists !!!!!\n\n";exit 1;}

if ($ResultPath=~/(.+)\/$/) { $ResultPath = "$1";}

`if [ ! -d "$ResultPath" ] ;then mkdir $ResultPath ; fi`;

sub display_help { usage(0); }





###############################
# 
# read vcf position
# 
open (INvcf,"$vcf_file")||die "$!";

my @tmp;
my ($pos,$ref,$alt);
my ($ref_len,$alt_len);
my %vcf;
my $vcf;
my $vcf_value;
my $snp;
my %snp;
my $insertion;
my %insertion;
my $deletion;
my %deletion;
my %type;
my $chrmosome;
my %chrmosome;
my %protein_coding;

while(<INvcf>){
	
	if (/^(.+?)\t/){
		
		#print $1;
		@tmp= split "\t",$_;
		
		$chrmosome=$tmp[0];
		$pos=$tmp[1];
		$ref=$tmp[3];
		$alt=$tmp[4];
		$chrmosome{$pos}=$chrmosome;
		
		$ref_len=length($ref);
		$alt_len=length($alt);
		
		$vcf_value=join "\t",$ref,$alt;
		$vcf{$pos}=$vcf_value;
		$protein_coding{$pos}=0;
		#print"$pos\t$vcf{$pos}\n";
		
		
		$vcf++;
		if ($ref_len == $alt_len){
			
			$snp{$pos} = $vcf_value;
			$type{$pos}="SNP";
			$snp++;
			
		}elsif($ref_len < $alt_len){
			
			$insertion{$pos}=$vcf_value;
			$type{$pos}="insertion";
			$insertion++;
			
		}elsif($ref_len > $alt_len){
			
			$deletion{$pos} = $vcf_value;
			$type{$pos}="deletion";
			$deletion++;
		}
	}
}
close INvcf;



###############################
# 
# genomic cds position
# 
open (INcds,"$gff")||die "$!";

my %locus_tag;		# record cds position
my $locus;
my $position;
my $position1;
my $position2;
my %protein_ID;		# record locus => protein ID & protein name
my $protein_ID;
my $protein_value;
my $i;
my $k;
my $j;
my $count=0;
my $inter=0;
my $pos_tmp;
my @pos_tmp=keys %vcf;
my @cds_locus;
my %output_format;
my %locus_tag_chrmosome;
my %gene_id;
my $gene_id;
my %non_coding_gene;
my $non_coding_gene;

###################
## get cds position
while(<INcds>)
{
	chomp;
	
	if (/^\#/){}
	else{
		@tmp= split "\t",$_;
		
		$chrmosome=$tmp[0];
		
		if($tmp[2]=~/gene/ )
		{
			
			if ($tmp[8]=~/\;old_locus_tag=(.+?)\;/)
			{
				$locus = $1;

			}
			elsif($tmp[8]=~/\;old_locus_tag=(.+)/)
			{
				$locus = $1;
			}
			
			#push @cds_locus,$locus;
			$locus_tag_chrmosome{$locus}=$chrmosome;
			#print "$locus\n";
			
			
			if ($tmp[8]=~/\;Dbxref=(.+?)\;/)
			{
				$gene_id{$1}=$locus;
			}
			
			##############################
			# read cds postion
			$position="$tmp[3]-$tmp[4]";
			$locus_tag{$locus}=$position;
			
			##############################
			# read cds info
			if(/\;gene_biotype=pseudogene\;/)
			{
				$protein_value="\tpseudogene";
				$protein_ID{$locus}=$protein_value;
				$non_coding_gene{$locus}="yes";
			}

		}
		elsif($tmp[2]=~/CDS/ || $tmp[2]=~/RNA/ )
		{
			#################
			#
			# get CDS info
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
			
			if ( $tmp[2]=~/RNA/ )
			{
				$non_coding_gene{$locus}="no";
			}
			else
			{
				$non_coding_gene{$locus}="yes";
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






##########################################################
# 
# Checking variants position 
# 
my @all_chromosomes=values %chrmosome;
my %cds_list;
my $variants_pos;
my %output_result;

foreach $locus (sort keys %locus_tag)
{
	## chromosome name:	$locus_tag_chrmosome{$locus}
	
	if ($locus_tag{$locus} =~ /([0-9]+)-([0-9]+)/)
	{
		$position1=$1;
		$position2=$2;
	}
	
	foreach $variants_pos (sort keys %vcf)
	{
		## chromosome name:	$chrmosome{$variants_pos}
		
		if ($locus_tag_chrmosome{$locus}==$chrmosome{$variants_pos} )
		{
			if ($position1 <= $variants_pos && $variants_pos <= $position2  )
			{
				$output_result{$variants_pos}="$variants_pos\t$type{$variants_pos}\t$vcf{$variants_pos}\t$locus\t$locus_tag{$locus}\t$non_coding_gene{$locus}\t$protein_ID{$locus}\n";
				$protein_coding{$variants_pos}++;
			}
			elsif ($position1 >= $variants_pos && $variants_pos >= $position2)
			{	
				$output_result{$variants_pos}="$variants_pos\t$type{$variants_pos}\t$vcf{$variants_pos}\t$locus\t$locus_tag{$locus}\t$non_coding_gene{$locus}\t$protein_ID{$locus}\n";
				$protein_coding{$variants_pos}++;
			}
		}
	}
}



#########################
# 
# Ouputing result form
# 
# 

open (OUT,">$ResultPath")||die "$!";

print OUT "Chomosome\t";
print OUT "variation position\t";
print OUT "variation type\t";
print OUT "REF\t";
print OUT "ALT\t";
print OUT "cds_locus_tag\t";
print OUT "cds_start-end\t";
print OUT "translation\t";
print OUT "protein_ID\t";
print OUT "protein_name\n";


my $noncoding=0;
foreach (sort keys %protein_coding)
{
	if( $protein_coding{$_} >0 )
	{
		print OUT "$chrmosome{$_}\t$output_result{$_}";
		
	}
	elsif($protein_coding{$_}==0)
	{
		print OUT "$chrmosome{$_}\t$_\t$type{$_}\t$vcf{$_}\tnon-coding region region\n";
		$noncoding++;
	}
	else{}	
}
close OUT;

my @array=keys %vcf;
my $array=@array;
print "CDS variant report:\n";
print "total variants : $array\n";
print "SNP : $snp\n";
print "insertion : $insertion\n";
print "deletion : $deletion\n";
print "noncoding region : $noncoding\n\n";

