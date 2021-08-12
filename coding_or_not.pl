#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
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
This program checking variant in coding or non-coding region:

Usage: $PROG [options] 
Options:
Input/Output:
  
  --input,-i STR          Input snp position annotation file [ default: variantion_all.txt ]
  
  --ResultPath,-o  STR        output file 
  					[ coding output default:  	  mutation_in_coding_region.txt     ]
					[ non-coding output default:  mutation_in_non-coding_region.txt ]
  
  --help,-h               Print this message
  --version,-v            Print version information

example:

perl $PROG \
-i variantion_all.txt \
-o ./result/test_sample \

EOF
    exit $exit_code;
}

GetOptions(
    "h" => \&display_help,
    "v" => \&display_version,
    "i=s" => \$input,
    "o=s" => \$ResultPath,
    
    "help" => \&display_help,
    "version" => \&display_version,
    "input=s" => \$input,
    "ResultPath=s" => \$ResultPath,
);
# Check input exists
if ( ! defined $input ) { print STDERR "\n !!!!! please input variantion_all.txt !!!!!! \n"; usage(0);}
if ( ! defined $ResultPath ) { $ResultPath = "variantion_all.txt";}


## Check input file exists
if ( ! -f "$input" ) { print STDERR "\n\n!!!! $input : variantion_all.txt file is not exists !!!!!\n\n";exit 1;}

if ($ResultPath=~/(.+)\/$/) { $ResultPath = "$1";}

`if [ ! -d "$ResultPath" ] ;then mkdir $ResultPath ; fi`;

sub display_help { usage(0); }


my $tmp;
my @tmp;
my $cds_name;
my %coding;
my %noncoding;
my $header;
my $chrmosome;
my $mutation_pos;

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
	
	$chrmosome=$tmp[0];

	if ($tmp[1]==""){}
	elsif($tmp[0]=~/^#/){}
	elsif ($tmp[5]=~/non-coding region/)
	{
		
		$noncoding{"$tmp[0]\t$tmp[1]"}=$_;
		
	}
	else
	{
		if ($chrmosome eq $tmp[0] )
		{
			if ($tmp[2] eq "SNP" )
			{
				$mutation_pos=$tmp[1];
				$coding{"$tmp[0]\t$tmp[1]"}=$_;
			}
			elsif( ( $tmp[1] - $mutation_pos ) < 10 )
			{
				#$coding{"$tmp[0]\t$tmp[1]"}=$_;
			}
			else
			{
				$mutation_pos=$tmp[1];
				$coding{"$tmp[0]\t$tmp[1]"}=$_;
			}
		}
		else
		{
			$mutation_pos=$tmp[1];
		}
	}

	

}


close IN;

open (OUTcoding,">$ResultPath/mutation_in_coding_region.txt")||die"$!";
print OUTcoding "$header";
foreach (sort keys %coding)
{
	print OUTcoding "$coding{$_}\n";
	
}

close OUTcoding;

open (OUTnoncoding,">$ResultPath/mutation_in_non-coding_region.txt")||die"$!";
print OUTnoncoding "$header";
foreach (sort keys %noncoding)
{
	print OUTnoncoding "$noncoding{$_}\n";
	
}

close OUTnoncoding;
