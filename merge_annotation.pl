#!/usr/bin/perl
use strict;
#############################################
# 
# This perl script merged mutation cds with uniprot annotation
# 
# usage : perl merge_annotation.pl $ARGV[0] $ARGV[1] > output.txt
# 
# $ARGV[0] >>> mutation cds file
# $ARGV[1] >>> uniprot annotation file
# 
# 
# 
# perl merge_annotation.pl HOL1_variantion_all.txt uniprot.txt > cds_merged.txt
# 
# uniprot.txt:
# -------------------------------------------------
# Entry	locus id	Entry name	Status	Protein names	Cross-reference (KO)	Gene ontology (biological process)	Gene ontology (cellular component)	Gene ontology (molecular function)
# A0A0M4D7G5	SPRI_0920	A0A0M4D7G5_STRPR	unreviewed	Heme chaperone HemW		porphyrin-containing compound biosynthetic process [GO:0006779]	cytoplasm [GO:0005737]	4 iron, 4 sulfur cluster binding [GO:0051539]; coproporphyrinogen oxidase activity [GO:0004109]; metal ion binding [GO:0046872]
# A0A0M4D695	SPRI_0952	A0A0M4D695_STRPR	unreviewed	Sugar-binding protein				carbohydrate binding [GO:0030246]; catalytic activity [GO:0003824]
# A0A0M4DNE1	SPRI_0989	A0A0M4DNE1_STRPR	unreviewed	Glycosyl transferase				transferase activity [GO:0016740]
# A0A0M4D1K4	SPRI_1048	A0A0M4D1K4_STRPR	unreviewed	ABC transporter ATP-binding protein				ATPase activity [GO:0016887]; ATP binding [GO:0005524]
# ------------------------------------------------
# 
# 
# 
# 
# Edit time:
# 20190908
# 
# 

my $tmp;
my @tmp;
my $cds_name;
my %cds;
my $header;

##############################
# 
# Doing uniprot annotation
# 
open (IN,"$ARGV[1]")||die "$!";
$header=<IN>;
while(<IN>)
{
	chomp;
	@tmp= split "\t",$_;
	$cds_name=$tmp[1];
	$cds{$cds_name}="$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]";
	

}


close IN;



##############################
# 
# Doing mutation cds file
# 
open (IN,"$ARGV[0]")||die "$!";
$header=<IN>;
print "Chomosome\t";
print "variation position\t";
print "variation type\t";
print "REF\t";
print "ALT\t";
print "cds_locus_tag\t";
print "cds_start-end\t";
print "translation\t";
print "ID\t";
print "product name\t";
print "product name(uniprot)\t";
print "KEGG (KO)\t";
print "Gene ontology (biological process)\t";
print "Gene ontology (cellular component)\t";
print "Gene ontology (molecular function)\n";

while(<IN>)
{
	chomp;
	@tmp= split "\t",$_;
	
	print "$_\t";
	foreach $cds_name (sort keys %cds)
	{
		if($tmp[5]=~/$cds_name/)
		{
			print "$cds{$cds_name}";
			
			}
	}
	print "\n";
	
}


close IN;











