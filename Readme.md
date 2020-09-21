# Genome mutation analysis

## This repo is doing:
###### Mutation analysis after `variant calling`
###### Extraction all mutation status from `vcf` file
###### **Check `amino acid mutation`**
###### **Check `non-coding region mutation`**
---


![image](https://github.com/HuangDeneil/genome-mutation-analysis/blob/master/image/workflow.png)


### Genome mutation analysis process after vcf file production step

### before vcf file production steps:

---

1) Qualuty control (`trimomatic`)

2) Reference genome alignmnet (`bwa`)

3) Alignment file convertion (optional) (`samtools`)

4) Output `vcf` file (`bcftools`)

`ps : vcf (variant calling file)`

<br>


#### After vcf file generation:

---


1. check_snp_position : This perl is checking mutation position with "gff file" from NCBI

``` bash
perl check_snp_position.pl sample.gff sample.vcf 
```
<br>

2. merge_annotation (need to formating `uniport` annotation file before this step)

``` bash
perl merge_annotation.pl HOL1_variantion_all.txt uniprot_annotation_info.txt > cds_merged.txt 
```
<br>

3. coding_or_not : This perl is `split` `coding` sequece mutation and `non-coding` sequece mutation

``` bash
perl coding_or_not.pl cds_merged.txt
```
<br>

4. noncoding_analysis : This perl is checking the nearest cds of non-coding mutation

``` bash
perl noncoding_analysis.pl sample.gff mutation_in_non-coding_region.txt
```
<br>

5. check_aa_sequence : This perl is primary check `amino acid mutation` with DNA sequence `translation` to amino acid sequence 

``` bash
perl check_aa_sequence.pl ATCC824_genome.fna codon_transfer.txt mutation_in_coding_region.txt 
```
<br>

6. aa_forword_check : 

``` bash
perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt 
```

<br>
##### More detail see Readme.sh

#### Comparing 3 mutation strain mutation genes:

```bash
7. 
perl cds_cross_check_v2.pl 0_S4_L001_final_report.txt 1_S5_L001_final_report.txt  2_S6_L001_final_report.txt
```



