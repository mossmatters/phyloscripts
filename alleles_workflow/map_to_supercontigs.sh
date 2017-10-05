#!/bin/bash

#PBS -q default
#PBS -j oe
#PBS -o ambiguity.out
#PBS -l nodes=1:ppn=1
#PBS -t 1-24

#cd $TMPDIR #Projects/artocarpus/alleles_paper/iupac_sequences
#prefix=$(tail -n $PBS_ARRAYID /home/mjohnson/Projects/artocarpus/alleles_paper/namelist_ajb.txt | head -1)

# This workflow will take the supercontig output of HybPiper and return a supercontig that 
# contains heterozygous positions as ambiguity bases. Uses paired reads.

#The script should be run on a FASTA file containing all the supercontigs of interest.


if [[ $# -eq 0 ]] ; then
    echo 'usage: hybpiper_ambiguity.sh supercontig.fasta readfile1.fq readfile2.fq'
    exit 1
fi

#########CHANGE THESE PATHS AS NEEDED###########

gatkpath=/opt/Software/GenomeAnalysisTK.jar
picardpath=/opt/Software/picard/build/libs/picard.jar

#############COMMAND LINE ARGUMENTS############

prefix=$1
read1fq=$2
read2fq=$3

mkdir $prefix
cd $prefix

#while read i
#do
#cat ~/Projects/artocarpus/alleles_paper/hybpiper/$prefix/$i/$prefix/sequences/intron/"$i"_supercontig.fasta
#done < ~/Projects/artocarpus/alleles_paper/newtargets_genelist.txt >> $prefix.supercontigs.fasta

supercontig=$prefix.supercontigs.fasta

#read1fq=~/Projects/artocarpus/alleles_paper/reads/"$prefix".R1.paired.fastq
#read2fq=~/Projects/artocarpus/alleles_paper/reads/"$prefix".R2.paired.fastq

#####STEP ZERO: Make Reference Databases

java -jar $picardpath CreateSequenceDictionary \
R=$supercontig 
bwa index $supercontig
samtools faidx $supercontig

#####STEP ONE: Map reads

echo "Mapping Reads"

bwa mem $supercontig $read1fq $read2fq | samtools view -bS - | samtools sort - -o $supercontig.sorted.bam

java -jar $picardpath FastqToSam  \
F1=$read1fq \
F2=$read2fq \
O=$supercontig.unmapped.bam \
SM=$supercontig

java -jar $picardpath MergeBamAlignment \
ALIGNED=$supercontig.sorted.bam \
UNMAPPED=$supercontig.unmapped.bam \
O=$supercontig.merged.bam \
R=$supercontig

#####STEP TWO: Mark duplicates

echo "Marking Duplicates"
java -jar $picardpath MarkDuplicates \
I=$supercontig.merged.bam \
O=$supercontig.marked.bam \
M=$supercontig.metrics.txt

#######STEP THREE: Identify variants, select only SNPs

echo "Identifying variants"

samtools index $supercontig.marked.bam
#samtools mpileup -B -f $supercontig $supercontig.marked.bam -v -u > $supercontig.vcf

java -jar $gatkpath \
-R $supercontig \
-T HaplotypeCaller \
-I $supercontig.marked.bam \
-o $supercontig.vcf



time java -jar $gatkpath \
-T SelectVariants \
-R $supercontig \
-V $supercontig.vcf \
-selectType SNP \
-o $supercontig.snps.vcf 


######STEP FOUR: Output new supercontig FASTA with ambiguity codes

echo "Generating IUPAC FASTA file"

java -jar $gatkpath \
-T FastaAlternateReferenceMaker \
-R $supercontig \
-o $supercontig.iupac \
-V $supercontig.snps.vcf \
-IUPAC $supercontig

cd ..
cp -r $prefix /home/mjohnson/Projects/artocarpus/alleles_paper/iupac_sequences/$prefix






