#!/usr/bin/perl
use strict;
##############################################################
# 
# Reorganized fasta file of sequence next line 
# 
# usage : perl fasta_reorganized.pl $ARGV[0] $ARGV[1] 
# 
# perl fasta_reorganized.pl ATCC824.fna ATCC824_reorganized.fna 
# 
# 
# 	$ARGV[0] >>> reference fasta (ATCC824.fna)
# 	$ARGV[1] >>> output file name
# 
# 
# edit time :
# 20190913
# 
# 


 my $input2;
 my $count2;
 my $seq2;
 my $j;
 my $sum;
 my @list;
 my @tmp2;
 my @namedir=qw//;
 my @tmpseq;
 my @seq;

 $input2=$ARGV[0] or die;
 
open (in2, "<$input2") || die"can't open  :$!";

print "working.....";
$count2=0;
 while(<in2>){
			chomp;
	 if($_=~/>/){
			@tmp2=split "\s",$_;
			$_=$tmp2[0];
			push @namedir,$_;

#eslf ($_=~/len/)   eslf ($_=~/path/)   eslf ($_=~/[(n*:n*-n*)]/) 

			 if($count2!=0){
					$seq2=join '',@tmpseq;
					push @seq,$seq2;
					@tmpseq=();
				 }
	 $count2++;
	 }else{
				push @tmpseq,$_;
				}
	}

		$seq2=join '',@tmpseq;
		push @seq,$seq2;
		@tmpseq=();
		print "stick last one \n";

open out, ">$ARGV[1]" || die"can't open  :$!";
 print "sort seq to file...\n";
 $sum=@namedir;
 print "\n\n>>you got $sum!\n";
 print "\n >>loading into result.fa .....\n\n";
 
for($j=0; $j<=@namedir;$j++)
{
	printf out "$namedir[$j]\t$seq[$j]\n"; 
}

close out;
