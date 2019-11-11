HOL1 genome mutation analysis

1. bwa
2. samtools
3. bcftools output vcf file
4. perl check_snp_position.pl GCF_000008765.1_ASM876v1_genomic.gff HOL_S10_L001.bwa.sorted.calls.filted_005.vcf 
5. perl merge_annotation.pl HOL1_variantion_all.txt uniprot.txt > cds_merged.txt
6. perl coding_or_not.pl cds_merged.txt
7. perl noncoding_analysis.pl GCF_000008765.1_ASM876v1_genomic.gff mutation_in_non-coding_region.txt
8. perl check_aa_sequence.pl ATCC824_reorganized.fna codon_transfer.txt mutation_in_coding_region.txt 
9. perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt



perl check_snp_position.pl GCF_000008765.1_ASM876v1_genomic.gff HOL_S10_L001.bwa.sorted.calls.filted_005.vcf 
perl merge_annotation.pl HOL1_variantion_all.txt uniprot.txt > cds_merged.txt
perl coding_or_not.pl cds_merged.txt
perl noncoding_analysis.pl GCF_000008765.1_ASM876v1_genomic.gff mutation_in_non-coding_region.txt
perl check_aa_sequence.pl ATCC824_reorganized.fna codon_transfer.txt mutation_in_coding_region.txt 
perl aa_forword_check.pl amino_acid_primary_check.txt mutation_cds.faa mutation_cds_modified.faa mutation_in_coding_region.txt



