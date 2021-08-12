#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $now;
my $PROG = basename $0;
my $exit_code;
my ($input, $input1, $input2, $input3, $input4, $output, $output1, $output2, $output3, $output4);
my ($input_path, $input_path1, $input_path2, $index_path, $ResultPath, $outpath, $thread);
my ($gtf_path, $genome_fna, $path_trimmomatic);
my $reads_status;

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
    my $def_thread_ct = exists $ENV{"NUM_THREADS"} ? (0 + $ENV{"NUM_THREADS"}) : 2;
    print STDERR <<EOF;
This program is doing variant calling, including three steps:
1. Trimmoatic           >>> reads quality control
2. bwa                  >>> genome alignment
3. samtools / bcftools  >>> alignment result sort/variant calling

Usage: $PROG [options] 
Options:
Input/Output:
  
  --fq1  STR              Path of input fastq folder path (input gzip file)
  --fq2  STR              Path of input fastq folder path (input gzip file)

  --trimpath STR          Path of trimmomatic
  --dbpath,-db  STR       bwa-mem2 index path
  --genome,-g STR         Reference genome fasta
  
  --TmpPath,-tmp  STR     output tmp file path [ default: ./tmp ]
  --ResultPath,-o  STR    output result path [ default: ./result ]
  --threads,-t NUM        Number of threads [default: $def_thread_ct ]
  
  --help,-h               Print this message
  --version,-v            Print version information

example:

perl do_variant_calling.pl \
-fq1 HOL_S10_L001_R1_001.fastq \
-fq2 HOL_S10_L001_R2_001.fastq \
-db ./library/db/bacteria \
-g ./library/genome.fna \
-tmp ./tmp/test_sample \
-o ./result/test_sample \
-t 16
EOF
    exit $exit_code;
}

GetOptions(
    "h" => \&display_help,
    "v" => \&display_version,
    "fq1=s" => \$input_path1,
    "fq2=s" => \$input_path2,
    "db=s" => \$index_path,
    "g=s" => \$genome_fna,
    "tmp=s" => \$outpath,
    "o=s" => \$ResultPath,
    "t=i" => \$thread,
    
    "help" => \&display_help,
    "version" => \&display_version,
    "trimpath=s" => \$path_trimmomatic,
    "dbpath=s" => \$index_path,
    "genome=s" => \$genome_fna,
    "TmpPath=s" => \$outpath,
    "ResultPath=s" => \$ResultPath,
    "threads=i" => \$thread,
);
# Check input exists
if ( ! defined $thread ) {$thread = $ENV{"NUM_THREADS"} || 2;}
if ( ! defined $input_path1 ) { print STDERR "\n!!!!! No reads input !!!!!! \n"; usage(0); }
if ( ! defined $input_path2 ) { $input_path2 = ""; }
if ( ! defined $outpath ) { $outpath = "tmp";}
if ( ! defined $ResultPath ) { $ResultPath = "result";}
if ( ! defined $index_path ) { print STDERR "\n !!!!! Index path is empty !!!!!! \n"; usage(0);}
if ( ! defined $genome_fna ) { print STDERR "\n !!!!! please input reference genome fasta !!!!!! \n"; usage(0);}
if ( ! defined $path_trimmomatic ) { print STDERR "\n\n !!!!! please define trimmomatic path !!!!!! \n"; usage(0);}


## Check input file exists
if ( ! -e "${index_path}.bwt.2bit.64" ) { print STDERR "\n\t--dbpath ${index_path}\n\n\t!!!! Index file is not exist !!!!\n\n"; exit 1;}
if ( ! -f "$path_trimmomatic" ) { print STDERR "\n\n\t!!!!! Cannot find program at $path_trimmomatic '!!!!\n\n"; exit 1;}
if ( ! -f "$genome_fna" ) { print STDERR "\n\n!!!! $genome_fna : genome fasta file is not exists !!!!!\n\n";exit 1;}
if ( ! -f "$input_path1" ) { print STDERR "\n\n!!!! $input_path1 : fastq folder is not exists !!!!!\n\n";exit 1;}
if ( ! $input_path2 eq "" ) { if( ! -f "$input_path2") { print STDERR "\n\n!!!! $input_path2 : fastq folder is not exists !!!!!\n\n";exit 1;}  }

if ( ! $input_path2 eq "" ){$reads_status = "paired-end";}
else{$reads_status = "single-end";}


# if ($input_path1=~/(.+)\/$/) { $input_path1 = "$1";}
# if ($input_path2=~/(.+)\/$/) { $input_path2 = "$1";}
if ($outpath=~/(.+)\/$/) { $outpath = "$1";}
if ($ResultPath=~/(.+)\/$/) { $ResultPath = "$1";}
if ($path_trimmomatic=~/(.+)\/$/) { $path_trimmomatic = "$1";}


`if [ ! -d "$outpath" ] ;then mkdir $outpath ; fi`;
`if [ ! -d "$ResultPath" ] ;then mkdir $ResultPath ; fi`;

sub display_help { usage(0); }


my $opt;
my ($i, $k, @tmp);
my ($cmd, $cmd1, $cmd2);
my (@file, $file, $file1, $file2);
my $index_name;

@tmp = split "/", ${index_path};
$index_name = $tmp[-1];


### Check fq file format (fq, fq.gz, fastq, fastq.gz)
my $type = "fq";

# $input_path =$input_path1;
if ($input_path1=~/(.+)\/(.+)$/) { $input_path = "$1"; $file = "$2";}
else{$input_path = "./"; $file = "$input_path1";}

@file=`cd $input_path ; ls *.fq  2>/dev/null`;
$i = @file;
if ($i == 0)
{
    @file = `cd $input_path ; ls *.fq.gz  2>/dev/null`;
    $type = "fq.gz";
    $i = @file;
    if ($i == 0)
    {
        @file = `cd $input_path ; ls *.fastq.gz  2>/dev/null`;
        $type = "fastq.gz";
        $i = @file;
        if ($i == 0)
        {
            @file = `cd $input_path ; ls *.fastq  2>/dev/null`;
            $type = "fastq";
            if ($i == 0){print STDERR "\n!!!input format error!!!\n";  usage(0); }
        }
    }
}

if ($type eq "fq" ){$file =~s/.fq//;}
if ($type eq "fq.gz" ){$file =~s/.fq.gz//;}
if ($type eq "fastq.gz" ){$file =~s/.fastq.gz//;}
if ($type eq "fastq" ){$file =~s/.fastq//;}



if ($reads_status eq "paired-end" )
{
    $file1 = "$input_path1";
    $file2 = "$input_path2";
    
    
    ## trimmomatic
    $input1="$outpath/$file1";
    $input2="$outpath/$file2";
    $output1="$outpath/$file.1.qtrim.fq";
    $output2="$outpath/$file.1.unparied.fq";
    $output3="$outpath/$file.2.qtrim.fq";
    $output4="$outpath/$file.2.unparied.fq";
    $opt="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
    $cmd = "java -jar $path_trimmomatic PE -threads $thread -phred33 $input1 $input2 $output1 $output2 $output3 $output4 $opt";


    ## bwa
    $input1 = "$outpath/$file1.qtrim.fq";
    $input2 = "$outpath/$file2.qtrim.fq";
    $output = "$outpath/$file.bwa.$index_name.sam";
    $cmd = "bwa mem $index_path $input1 $input2 -o $output";


    ## samtools
    ## sam to bam
    $input = "$outpath/$file.bwa.$index_name.sam";
    $output = "$outpath/$file.bwa.$index_name.sort.bam";
    $cmd = "samtools sort -@ $thread $input | samtools view -@ $thread -b -F 4  > $output";

    ##samtools index  ( 產生 .bam.bai )
    $input = "$outpath/$file.bwa.$index_name.sort.bam";
    $cmd = "samtools index -b $input";

    #### creating VCF/BCF file
    $input = "$outpath/$file.bwa.$index_name.sort.bam";
    $output = "$outpath/$file.bwa.$index_name.calls.bcf";
    $cmd = "bcftools mpileup -f ./index/reference.fna $input | bcftools call -mv -Ob -o calls.bcf";

    ##quality score (the QUAL) column, Q=-10*log(P value)
    $input = "$outpath/$file.bwa.$index_name.calls.bcf";
    $output = "$outpath/$file.bwa.$index_name.filted.vcf";
    $cmd = "bcftools view -i '\%QUAL>=20' $input -Ov -o $output";
}
elsif ($reads_status eq "single-end" )
{
    
    $file1 = "$input_path1";
        
        
    ## trimmomatic
    $input1="$outpath/$file1.fq";
    $output1="$outpath/$file.qtrim.fq";
    $output2="$outpath/$file.unparied.fq";
    $opt="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
    $cmd = "java -jar $path_trimmomatic SE -threads $thread -phred33 $input1 $output1 $output2 $opt";


    ## bwa
    $input1 = "$outpath/$file.qtrim.fq";
    $output = "$outpath/$file.bwa.$index_name.sam";
    $cmd = "bwa mem $index_path $input1 -o $output";


    ## samtools
    ## sam to bam
    $input = "$outpath/$file.bwa.$index_name.sam";
    $output = "$outpath/$file.bwa.$index_name.sort.bam";
    $cmd = "samtools sort -@ $thread $input | samtools view -@ $thread -b -F 4  > $output";

    ##samtools index  ( 產生 .bam.bai )
    $input = "$outpath/$file.bwa.$index_name.sort.bam";
    $cmd = "samtools index -b $input";

    #### creating VCF/BCF file
    $input = "$outpath/$file.bwa.$index_name.sort.bam";
    $output = "$outpath/$file.bwa.$index_name.calls.bcf";
    $cmd = "bcftools mpileup -f $genome_fna $input | bcftools call -mv -Ob -o calls.bcf";

    ##quality score (the QUAL) column, Q=-10*log(P value)
    $input = "$outpath/$file.bwa.$index_name.calls.bcf";
    $output = "$outpath/$file.bwa.$index_name.filted.vcf";
    $cmd = "bcftools view -i '\%QUAL>=20' $input -Ov -o $output";
}

