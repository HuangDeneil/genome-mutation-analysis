# Genome mutation analysis

## This repo is doing mutation analysis
---


[![image](https://github.com/HuangDeneil/genome-mutation-analysis/blob/master/image/workflow.png)]


### Genome mutation analysis process after vcf file production step

### before vcf file production steps:
##### 1.Alignmnet (bwa)
##### 2.Alignment file convertion (optional) (samtools)
##### 3.Output vcf file (bcftools)


<h4>After vcf file created:</h4>
<p>1. perl check_snp_position.pl sample.gff sample.vcf 
<p>2. perl merge_annotation.pl HOL1_variantion_all.txt uniprot_annotation_info.txt > cds_merged.txt </p>
<p>3. perl coding_or_not.pl cds_merged.txt</p>
<p>4. perl noncoding_analysis.pl sample.gff mutation_in_non-coding_region.txt </p>
<p>5. perl check_aa_sequence.pl ATCC824_genome.fna codon_transfer.txt mutation_in_coding_region.txt </p>
<p>6. perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt </p>

<h2>more detail see Readme.sh</h2>

<h3>Comparing 3 mutation strain mutation genes:</h3>
<p>7. perl cds_cross_check_v2.pl 0_S4_L001_final_report.txt 1_S5_L001_final_report.txt  2_S6_L001_final_report.txt</br>




