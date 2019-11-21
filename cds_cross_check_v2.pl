#!/usr/bin/perl
use strict;

#############################################
# 
# This perl script sorting geranal table
# 
# usage : perl cds_cross_check.pl $ARGV[0] $ARGV[1] $ARGV[2] 
# 
# $ARGV[0] >>> 0_final_report.txt
# $ARGV[1] >>> 1_final_report.txt
# $ARGV[2] >>> 2_final_report.txt
# 
# perl cds_cross_check.pl 0_final_report.txt 1_final_report.txt 2_final_report.txt
# 
# perl cds_cross_check_v2.pl 0_S4_L001_final_report.txt 1_S5_L001_final_report.txt  2_S6_L001_final_report.txt
# 
# 
# 0_final_report.txt:
# --------------------------------------------------
# cds_locus_tag	cds_start-end	mutation_point	translation	aa_change_info	ID	product name	product name(uniprot)	EC	KEGG (KO)	Gene ontology (biological process)	Gene ontology (cellular component)	Gene ontology (molecular function)
# SPRI_0343	442448-442846	1	yes	animo acid change but maybe work	WP_053556652.1	DUF1992 domain-containing protein	Molecular chaperone DnaJ					
# SPRI_0360	461513-462394	1	yes	animo acid change but maybe work	WP_037772923.1	LLM class flavin-dependent oxidoreductase	Oxidoreductase				oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen [GO:0016705]	
# SPRI_0390	498294-498758	1	yes	mutilple or lose stop codon	WP_053556657.1	methylmalonyl-CoA mutase	Methylmalonyl-CoA mutase	K01849;			cobalamin binding [GO:0031419]; isomerase activity [GO:0016853]; metal ion binding [GO:0046872]	
# SPRI_0454	572417-572926	1	yes	animo acid change but maybe work	WP_053557676.1	SRPBCC family protein	Polyketide cyclase					
# 
# 
# 
# 
# 
# 
# 


if (!$ARGV[0] || !$ARGV[1] || !$ARGV[2] )
{
	print "\n\n";
	print "!!!!!!!No list detect!!!!\n";
	print "Uasge:\n";
	print "\tperl cds_cross_check.pl \$ARGV[0] \$ARGV[1] \$ARGV[2] \n\n";
	print "\$ARGV[0] >>> 0_final_report.txt\n";
	print "\$ARGV[1] >>> 1_final_report.txt\n";
	print "\$ARGV[2] >>> 2_final_report.txt\n";
	die "error";
}

my @faafile=@ARGV;
my $faafile;
my $i;
my $file;
my $key;
my %data;
my %cross_check;
my @tmp;
my %aa_change_info;
my $count=0;
my %seq;
my $header;

for ( $i=0; $i<@ARGV; $i++ )
{
	###############################
	# 
	# Reading final report file
	# 
	$file = $ARGV[$i];
	open(IN,"$file")||die "$!";
	$header= <IN>;
	while (<IN>)
	{
		chomp;
		@tmp = split "\t",$_;
		if ($tmp[0] eq ""){}
		else{
			$key = "$file\t$tmp[0]";
			$data{$key} = $_;
			
			if ( $tmp[4] =~ /normal/ ){}
			else
			{
				$aa_change_info{$key}=$tmp[4];
				$cross_check{$tmp[0]}="$cross_check{$tmp[0]}\t$file";
			}
			
		}
	}
	close IN;
	
	
	$faafile=$ARGV[$i];
	$faafile =~s/_final_report.txt/_mutation_cds.faa/;
	###############################
	# 
	# Reading mutation amino acid file
	# 
	$file = $faafile;
	open(IN,"$file")||die "$!";
	
	while (<IN>)
	{
		chomp;
		
		if (/^>/)
		{
			@tmp = split "\t",$_;
			$tmp[0]=~s/>//;
			$key = "$ARGV[$i]\t$tmp[0]";
			$count++;
		}
		elsif( $count == 1 )
		{
			$seq{$key} = $_;
			$count=0;
		}
	}
	close IN;
	
}

my $str;
my $cds;
my ($a,$b,$c,$ab,$ac,$bc,$abc);
my (@a,@b,@c,@ab,@ac,@bc,@abc);
my (%a,%b,%c,%ab,%ac,%bc,%abc);
my ($key0,$key1,$key2,$key3,$keya,$keyb,$keyc,$keyabc);

foreach $cds(sort keys %cross_check)
{
	######
	#  0
	if ( $cross_check{$cds} eq "\t$ARGV[0]" )
	{
		$a{$cds}=1;
	}
	
	######
	#  1
	if ( $cross_check{$cds} eq "\t$ARGV[1]" )
	{
		$b{$cds}=1;
	}
	
	######
	#  2
	if ( $cross_check{$cds} eq "\t$ARGV[2]" )
	{
		$c{$cds}=1
	}
	
	#########
	#  0 & 1
	if ( $cross_check{$cds} eq "\t$ARGV[0]\t$ARGV[1]" )
	{
		$keya="$ARGV[0]\t$cds";
		$keyb="$ARGV[1]\t$cds";
		if ( $seq{$keya} eq $seq{$keyb} )  ## a & b aa seq same  >> ab unique
		{
			$ab{$cds}=1;	
			
		}
		else		## a & b aa seq different  >> a unique & b unique
		{
			$a{$cds}=1;
			$b{$cds}=1;
		}
		
	}
	
	#########
	#  0 & 2
	#$str=;
	if ( $cross_check{$cds} eq "\t$ARGV[0]\t$ARGV[2]" )
	{
		$keya="$ARGV[0]\t$cds";
		$keyc="$ARGV[2]\t$cds";
		#print "$cross_check{$cds}\n";
		
		if ( $seq{$keya} eq $seq{$keyc} )  ## a & c aa seq same >> ac unique
		{
			$ac{$cds}=1;
			
		}
		else		## a & c aa seq different>> a uniuq & c unique
		{
			$a{$cds}=1;
			$c{$cds}=1;
			
			
		}
	}
	
	
	#########
	#  1 & 2
	if ( $cross_check{$cds} eq "\t$ARGV[1]\t$ARGV[2]" )
	{
		$keyb="$ARGV[1]\t$cds";
		$keyc="$ARGV[2]\t$cds";
		if ( $seq{$keyb} eq $seq{$keyc} )  ## b & c aa seq same >>  bc unique
		{
			$bc{$cds}=1;
		}
		else		## b & c aa seq different >> b unique & c unique
		{
			$b{$cds}=1;
			$c{$cds}=1;
		}	
	}
	
	
	############
	#  0 & 1 & 2
	
	if ( $cross_check{$cds} eq "\t$ARGV[0]\t$ARGV[1]\t$ARGV[2]" )
	{
		$keya="$ARGV[0]\t$cds";
		$keyb="$ARGV[1]\t$cds";
		$keyc="$ARGV[2]\t$cds";
		
		#$abc{$cds}=1;
		if ( $seq{$keya} eq $seq{$keyb} )  ## a & b aa seq same
		{
			if ( $seq{$keya} eq $seq{$keyc} ) ## a & b & c aa seq same
			{
				$abc{$cds}=1;
			}
			else ## a & b aa seq same but c differnt
			{
				$ab{$cds}=1;
				$c{$cds}=1;
			}
		}
		elsif ( $seq{$keya} eq $seq{$keyc} )  ## a & c aa seq same
		{
			if ( $seq{$keyb} eq $seq{$keyc} ) ## a & b & c aa seq same
			{
				$abc{$cds}=1;
			}
			else ## a & c aa seq same but b differnt
			{
				$ac{$cds}=1;
				$b{$cds}=1;
			}
		}
		elsif ( $seq{$keyb} eq $seq{$keyc} )  ## b & c aa seq same
		{
			if ( $seq{$keya} eq $seq{$keyc} ) ## a & b & c aa seq same
			{
				$abc{$cds}=1;
			}
			else ## b & c aa seq same but a differnt
			{
				$bc{$cds}=1;
				$a{$cds}=1;
			}
		}
		else		## a & b aa seq different
		{
			$a{$cds}=1;
			$b{$cds}=1;
			$c{$cds}=1;
		}
			
	}
} 
@a=sort keys %a;
@b=sort keys %b;
@c=sort keys %c;
@ab=sort keys %ab;
@ac=sort keys %ac;
@bc=sort keys %bc;
@abc=sort keys %abc;

$a=@a;
$b=@b;
$c=@c;
$ab=@ab;
$ac=@ac;
$bc=@bc;
$abc=@abc;

print "0 ($a)\n";
print "1 ($b)\n";
print "2 ($c)\n";
print "01 ($ab)\n";
print "02 ($ac)\n";
print "12 ($bc)\n";
print "012 ($abc)\n";

open(OUT, ">a.txt")||die "$!";
print OUT "$header";
foreach (@a)
{
	$key="$ARGV[0]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">b.txt")||die "$!";
print OUT "$header";
foreach (@b)
{
	$key="$ARGV[1]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">c.txt")||die "$!";
print OUT "$header";
foreach (@c)
{
	$key="$ARGV[2]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">ab.txt")||die "$!";
print OUT "$header";
foreach (@ab)
{
	$key="$ARGV[0]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">ac.txt")||die "$!";
print OUT "$header";
foreach (@ac)
{
	$key="$ARGV[0]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">bc.txt")||die "$!";
print OUT "$header";
foreach (@bc)
{
	$key="$ARGV[1]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

open(OUT, ">abc.txt")||die "$!";
print OUT "$header";
foreach (@abc)
{
	$key="$ARGV[0]\t$_";
	print OUT "$data{$key}\n";
}
close OUT;

