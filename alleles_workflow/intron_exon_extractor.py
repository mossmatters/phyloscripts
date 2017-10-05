#!/usr/bin/env python
# Script to use the GFF files from intronerate and the ambiguity-encoded FASTA files to generate separate intron and exon files for each gene.

import re,sys,os,errno,shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

my_re = re.compile(r"([0-9]+ )(.+):1")

#Fix the names in IUPAC file

prefix=sys.argv[1]

if os.path.isdir(prefix):
    os.chdir(prefix)

if os.path.isfile("{}.supercontigs.fasta.iupac".format(prefix)):
    with open("{}.supercontigs.iupac.fasta".format(prefix),'w') as outfile:
        for seq in SeqIO.parse("{}.supercontigs.fasta.iupac".format(prefix),'fasta'):
            seq.id = my_re.sub("\g<2>",seq.description)
            seq.description = ''
            SeqIO.write(seq,outfile,'fasta')
        

#Parse GFF into dictionaries for each gene (one for introns, one for exons)
# ASSUMES THE GFF IS SORTED WITHIN EACH GENE!!!!

intron_dict = {}
exon_dict = {}

gff_fn = prefix+".intronerate.gff"#sys.argv[2]
for line in open(gff_fn):
    line=line.split()
    if line[2] == "exon":
        try:
            exon_dict[line[0]].append((int(line[3])-1,int(line[4])))
        except KeyError:
            exon_dict[line[0]] = [(int(line[3])-1,int(line[4]))]
#    elif line[2] == "intron":
#        try:
#            intron_dict[line[0]].append((int(line[3])-1,int(line[4])))
#        except KeyError:
#            intron_dict[line[0]] = [(int(line[3])-1,int(line[4]))]

try:
    supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.iupac.fasta".format(prefix),'fasta'))
    dataType = "iupac"
except IOError:
    try:
        supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.alleles.fasta".format(prefix),'fasta'))
        dataType = "alleles"
    except IOError:
        try:
            supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.svdq.fasta".format(prefix),'fasta'))
            dataType = "svdq"
        except IOError:
            supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.default.fasta".format(prefix),'fasta'))
            dataType = 'default'
        
    
for gene in exon_dict:
    try:
        geneLength = len(supercontig_dict[gene])
    except KeyError:
        haploGeneName = "{}_h1-{}".format(prefix,gene.split("-")[-1])
        geneLength = len(supercontig_dict[haploGeneName])
    exon_ranges = exon_dict[gene]
#    intron_dict[gene] = [(0,exon_dict[0][0]),exon_dict[0][1]]
    for exon_interval in range(len(exon_ranges)+1):
        if exon_interval == 0:
            intron_dict[gene] = [(0,exon_ranges[exon_interval][0]-1)]
        elif exon_interval == len(exon_ranges)  :
            intron_dict[gene].append((exon_ranges[-1][1],geneLength))

        else:
            start = exon_ranges[exon_interval - 1][1] 
            stop = exon_ranges[exon_interval][0] - 1
            intron_dict[gene].append((start,stop))
#print(intron_dict["NZ866-gene001.single"])
    
    

            

newseq = ''

for seqType in ["exon","intron","supercontig"]:
    if os.path.exists(seqType):
        shutil.rmtree(seqType)
    os.makedirs(seqType)

for gene in supercontig_dict:
    if gene.startswith("McBryde"):
        geneName = gene.split("-")[2]
        sampleName = gene.split("-")[1]
    else:
        geneName = gene.split("-")[1]
        sampleName = gene.split("-")[0]
            
    with open("exon/{}.{}.FNA".format(geneName,dataType),'a') as exonout:
        newseq = ''
        exonLookupName = supercontig_dict[gene].id.replace("_h1",'')
        exonLookupName = exonLookupName.replace("_h2",'')
        if exonLookupName not in exon_dict:
            continue
        for gff_interval in exon_dict[exonLookupName]:
            newseq += supercontig_dict[gene].seq[gff_interval[0]:gff_interval[1]]
        exonout.write(">{}\n{}\n".format(sampleName,newseq))
    with open("intron/{}.intron.{}.fasta".format(geneName,dataType),'a') as intronout:
        newseq=''
        for gff_interval in intron_dict[exonLookupName]:
            newseq += supercontig_dict[gene].seq[gff_interval[0]:gff_interval[1]]
        intronout.write(">{}\n{}\n".format(sampleName,newseq))
    
    with open("supercontig/{}.supercontig.{}.fasta".format(geneName,dataType),'a') as supercontigout:
        supercontigout.write(">{}\n{}\n".format(sampleName,supercontig_dict[gene].seq))
        