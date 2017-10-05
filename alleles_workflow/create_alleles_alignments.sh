#!/bin/bash

#PBS -q default
#PBS -l nodes=1:ppn=12
#PBS -j oe
#PBS -o iupac_alignments.out


# Shell script to recreate the IUPAC ambiguity coded alignments for Artocarpus.


cd ~/Projects/artocarpus/alleles_paper/haplotype_sequences

genelist=/home/mjohnson/Projects/artocarpus/alleles_paper/singlecopy_genelist.txt
namelist=/home/mjohnson/Projects/artocarpus/alleles_paper/namelist_ajb.txt

##########EXONS############

###### Exon sequences generated from HybPiper output:

mkdir -p exon
rm exon/* 
parallel "cat /home/mjohnson/Projects/artocarpus/alleles_paper/haplotype_sequences/{1}/exon/{2}.alleles.FNA >> exon/{2}.alleles.FNA" :::: $namelist :::: $genelist

##### Alignments with MACSE

parallel --eta macse -prog alignSequences -seq exon/{}.alleles.FNA :::: $genelist

##### Replace frame shifts ! with gaps -

mkdir -p macse
rm macse/*
mv exon/*_macse* macse

parallel sed -i "s/\!/-/g" macse/{}.alleles_macse_NT.fasta :::: $genelist

##### Trim alignments to retain only sites present in 75% of taxa

mkdir -p exon_trimmed
rm exon_trimmed/*

parallel "trimal -gt 0.75 -in macse/{}.alleles_macse_NT.fasta -out exon_trimmed/{}.alleles.macse.trimmed.FNA" :::: $genelist

#Fix McBryde-MV2

parallel sed -i -E 's/McBryde-MV2/McBryde/g' exon_trimmed/{}.macse.trimmed.FNA :::: $genelist

###########INTRONS#########

##### Intron sequences generated from HybPiper (intronerate.py):

mkdir -p intron
rm intron/*
parallel "cat /home/mjohnson/Projects/artocarpus/alleles_paper/haplotype_sequences/{1}/intron/{2}.intron.alleles.fasta >> intron/{2}.intron.alleles.fasta" :::: $namelist :::: $genelist

# Remove gene name from intron sequence files

#parallel sed -i -E 's/-.+$//g' intron/{}.intron.alleles.fasta :::: $genelist

# Align intron sequences with MAFFT. Timeout because of known huge sequence.
mkdir -p mafft
rm mafft/*

parallel --timeout 4000% --eta "mafft --maxiterate 1000 --globalpair --preservecase intron/{}.intron.alleles.fasta > mafft/{}.intron.alleles.mafft.fasta" :::: $genelist

# Trim alignments to retain only sites present in 75% of taxa

mkdir -p intron_trimmed
rm intron_trimmed/*

parallel "trimal -gt 0.75 -in mafft/{}.intron.alleles.mafft.fasta -out intron_trimmed/{}.intron.alleles.mafft.trimmed.fasta" :::: $genelist


# Combine alignments

#parallel python ../../combine_alignments.py exon_trimmed/{}.iupac.macse.trimmed.FNA intron_trimmed/{}.intron.iupac.mafft.trimmed.fasta {} :::: $genelist

#mkdir -p ../artocarpus_alignments/default/
#mv *.fasta ../artocarpus_alignments/default/
#mv *.partition ../artocarpus_alignments/default
