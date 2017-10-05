#!/bin/bash

#Script to prepare phased haplotype sequences for each for one sample. 

prefix=$1
genelist=genelist.txt
mkdir -p $prefix
cd $prefix
rm -r *

#Run bcftools to extract sequences

bgzip -c $prefix.supercontigs.fasta.snps.whatshap.vcf > $prefix.supercontigs.fasta.snps.whatshap.vcf.gz
tabix $prefix.supercontigs.fasta.snps.whatshap.vcf.gz
mkdir -p phased_bcftools
rm phased_bcftools/*

parallel "samtools faidx $iupac_dir/$prefix.supercontigs.fasta $prefix-{1} | bcftools consensus -H 1 $prefix.supercontigs.fasta.snps.whatshap.vcf.gz > phased_bcftools/$prefix-{1}.phased.fasta" :::: $genelist 
parallel "samtools faidx $iupac_dir/$prefix.supercontigs.fasta $prefix-{1} | bcftools consensus -H 2 $prefix.supercontigs.fasta.snps.whatshap.vcf.gz >> phased_bcftools/$prefix-{1}.phased.fasta" :::: $genelist 

cd ..


