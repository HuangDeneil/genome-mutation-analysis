#!/bin/bash
#HCCB10218 genome mutation analysis

########################################
# Uasge :
# 	bash Readme.sh $1
# 
# 	$1 >>> vcf file
# 
# bash Readme.sh 0_S4_L001.bwa.sorted.filtered_0.05.vcf

input_vcf_file="$1"


refernce_gff="./reference_info/GCF_000008765.1_ASM876v1_genomic.gff"
reference_genome_fasta="./reference_info/GCF_000008765.1_ASM876v1_genomic.fna"
ref_genome="ATCC824.fna"
codon_table="./codon_transfer.txt"

uniprot="uniprot.txt"

necessary_file="\
$input_vcf_file \
$refernce_gff \
$reference_genome_fasta \
$codon_table"

for i in $necessary_file
do
	if [ ! -f "$i" ] ; then
		echo "$i lose"
		exit 1
	fi
done

script_folder_path="./script"
code="\
check_snp_position.pl \
merge_annotation.pl \
coding_or_not.pl \
noncoding_analysis.pl \
fasta_reorganized.pl \
check_aa_sequence.pl \
aa_forword_check.pl"

for i in $code
do
	
	if [ ! -f "$i" ];then
		if [ ! -f "$script_folder_path/$i" ];then
			echo "$script_folder_path/$i lose"
			exit 1
		else
			ln -s $script_folder_path/$i	
		fi
	fi
done

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##################################
##								##
##		After vcf file 			##
##								##
##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


#################################
# 
# 1. primary vcf reoganized
# 
echo -e "primary vcf reoganization (check_snp_position.pl)\n\n"
perl check_snp_position.pl $refernce_gff $input_vcf_file  
# output file : variantion_all.txt


##############################
# 
# 2. reoganized cds mutation uniprot info
# 
echo -e "reoganized cds mutation uniprot info (merge_annotation.pl)\n\n"
if [ ! -f "variantion_all.txt" ] ;then 
	echo "Something wrong with check_snp_position.pl"
	exit 1
fi

if [ ! -f "$uniprot" ] ;then
	echo "Trying to download uniprot annotation file first" 
	echo "uniprot colomn need" 
	echo "locus_id	Entry	Entry name	Status	Protein names	Cross-reference (KEGG)	Gene ontology (biological process)	Gene ontology (cellular component)	Gene ontology (molecular function)	Gene ontology (GO)" 
# locus_id	
# Entry	
# Entry name	
# Status	
# Protein names	
# Cross-reference (KEGG)	
# Gene ontology (biological process)	
# Gene ontology (cellular component)	
# Gene ontology (molecular function)	
# Gene ontology (GO)

	exit 1
	
fi


perl merge_annotation.pl variantion_all.txt $uniprot > cds_merged.txt


##############################
# 
# 3. Splited coding or not into two file
# 
echo -e "Splited coding or not into two file (coding_or_not.pl)\n\n"
perl coding_or_not.pl cds_merged.txt 
# output file 1 : mutation_in_coding_region.txt 
# output file 2 : mutation_in_non-coding_region.txt
#


##############################
# 
# 4. Non-coding analysis 
# 
echo "Doing noncoding_analysis.pl"
if [ ! -f "mutation_in_non-coding_region.txt" ] ;then 
	echo -e "Something wrong with coding_or_not.pl\n\n"
	exit 1
fi
perl noncoding_analysis.pl $refernce_gff mutation_in_non-coding_region.txt


##############################
# 
# 5. Re-Organized geneome fasta file
# 
perl fasta_reorganized.pl $reference_genome_fasta $ref_genome 


##############################
# 
# 6. Modifying into mutated fasta file (DNA & amino acid sequence)
# 
echo -e "Doing check_aa_sequence.pl\n\n"

if [ ! -f "mutation_in_coding_region.txt" ] ;then 
	echo "Something wrong with coding_or_not.pl"
	exit 1
fi

perl check_aa_sequence.pl $ref_genome $codon_table mutation_in_coding_region.txt 
# Created folder : Sequence_comparison/
# inside folder:
# 	mutation_cds_modified.faa
# 	mutation_cds_modified.fna
# 	amino_acid_primary_check.txt


for i in $code
do
	unlink $i
done



##############################
# 
# 7. Comparison of sequece difference
# 
if [ -d "Sequence_comparison" ] ;then
	cd Sequence_comparison
	if [ ! -f "mutation_in_coding_region.txt" ] ;then
		ln -s ../mutation_in_coding_region.txt
	fi
fi

if [ ! -f "aa_forword_check.pl" ] ;then
	ln -s ../script/aa_forword_check.pl ./
fi


perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt
# ./Sequence_comparison/
# 	final_report.txt
# 	
mv final_report.txt ${input_vcf_file}_final_report.txt


