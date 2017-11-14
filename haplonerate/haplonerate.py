#!/usr/bin/env python


import sys,os,argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

helptext = '''This script takes two files, containing phased haplotype sequences for one
or more genes and edits the sequences to retain only variable sites in the largest phase
block. The two haplotype sequences can be generated by bcftools, for example.

There are three options for editing the sequences:

delete: retain only the longest phase block, delete the rest of the sequence
ref: use reference sequences to fill the rest of the sequence outside the longest block
mask: fill sequence not in the longest phase block with N

haplonerate reads the phase blocks from a GTF file, such as the one produced by the stats 
function in whatshap

'''

def get_gtf_dict(gtf_fn):
    gtf_dict = {}
    for line in open(gtf_fn):
        line = line.split()
        geneName = line[0]#.split("-")[-1]
        phase_range = (int(line[3]),int(line[4]))
        try:
            gtf_dict[geneName].append(phase_range)
        except KeyError:
            gtf_dict[geneName] = [phase_range]
    return gtf_dict



#prefix = sys.argv[1]

#Use the bcftools method to extract haplotypes from each gene

#seqdirectory = "/home/mjohnson/Projects/artocarpus/alleles_paper/iupac_sequences/{}".format(prefix)
#geneList = set([x.rstrip() for x in open("/home/mjohnson/Projects/artocarpus/alleles_paper/newtargets_genelist.txt")])
#os.chdir(prefix)

#Read in the ambiguity coded sequences into a dictionary

#iupac_dict = SeqIO.to_dict(SeqIO.parse("{}/{}.supercontigs.iupac.fasta".format(seqdirectory,prefix),'fasta'))



#Use the GTF from Whatshap to determine the longest phase block for the sequence.


def getLargestPhaseBlock(ranges,seqLength):
    '''Given the phase blocks for a sequence, return the most inclusive range'''
    longestblock = 0
    for r in range(len(ranges)):
        if r == 0:
            start = 1
            if len(ranges) > 1:
                end = ranges[r+1][0] - 1
            else:    
                end = seqLength
            
        elif r == len(ranges) - 1 :
            start = ranges[r-1][1] + 1
            end = seqLength
            
        else:
            start = ranges[r-1][1] + 1
            end = ranges[r+1][0] - 1
            
            
        if end - start > longestblock:
            most_inclusive_range = (start,end)
            longestblock = end - start
    #print seqLength,ranges,most_inclusive_range    
    return most_inclusive_range

def insertPhase(iupacSeq,haploSeq,phaseBlock,newSeqID):
    '''Given an IUPAC sequence, the phased haplotype sequence, and the longest phaseBlock,
    Return one sequence with phased characters in the block, IUPAC sequences outside it'''
    
    newSeq = ''
    for c in range(len(iupacSeq.seq)):
        if  phaseBlock[0] -1 <= c <= phaseBlock[1] - 1:
            newSeq += haploSeq.seq[c]
        else:
            newSeq += iupacSeq.seq[c]
    return SeqRecord(Seq(newSeq),id=newSeqID,description='')
    

def replace_with_ref(seq1,seq2,ref,phaseBlock):
    if seq1.seq == seq2.seq:
        return [SeqRecord(ref.seq,id=ref.id,description='')]
    else:
        haplo1 = insertPhase(ref,seq1,phaseBlock,"{}_h1".format(seq1.id))
        haplo2 = insertPhase(ref,seq2,phaseBlock,"{}_h2".format(seq2.id))
        return [haplo1,haplo2]

def delete_extra(seq1,seq2,ref,phaseBlock):
    if seq1.seq == seq2.seq:
        return [SeqRecord(ref.seq,id=ref.id,description='')]
    else:
        haplo1 = SeqRecord(seq1.seq[phaseBlock[0]:phaseBlock[1]],id="{}_h1".format(seq1.id),description='')
        haplo2 = SeqRecord(seq2.seq[phaseBlock[0]:phaseBlock[1]],id="{}_h2".format(seq2.id),description='')
        return [haplo1,haplo2]

def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("gtf",help="gtf file annotating the positions of phase blocks for each gene")
    parser.add_argument("haplotype_files",help="Two FASTA files containing sequences for one or more genes",nargs="+")
    parser.add_argument("--reference","-r",help="FASTA file of reference sequences, required with --edit ref")
    parser.add_argument("--edit",help="How to deal with sites outside longest phase block. Default: ref",default="ref",choices=["ref","delete","mask"])
    parser.add_argument("--output",'-o',help="Output FASTA containing haplotype sequences for each gene. default = stdout",default=sys.stdout)
    parser.add_argument("--block",help="file to write phase block information") 
    args = parser.parse_args()
    
    if len(args.haplotype_files) != 2:
        print("Please supply exactly two haplotype FASTA files!\n")
        sys.exit(1)
    
    gtf_dict = get_gtf_dict(args.gtf)
    
    #if args.edit == "ref":
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.reference,'fasta'))
    
    haplotype1_dict = SeqIO.to_dict(SeqIO.parse(args.haplotype_files[0],'fasta'))
    haplotype2_dict = SeqIO.to_dict(SeqIO.parse(args.haplotype_files[1],'fasta'))
    geneList = set(haplotype1_dict.keys())
    seqs_to_write = []
    phase_report = []
    for gene in geneList:
        if gene in gtf_dict:
            if gene in ref_dict:
                phaseBlock = getLargestPhaseBlock(gtf_dict[gene], len(ref_dict[gene]))
                phase_report.append("{}\t{}\t{}\t{}\t{}".format(gene,len(gtf_dict[gene]),len(ref_dict[gene]),phaseBlock[0],phaseBlock[1]))
                if args.edit == 'ref':
                    seqs_to_write += replace_with_ref(haplotype1_dict[gene],haplotype2_dict[gene],ref_dict[gene],phaseBlock)
                elif args.edit == "delete":
                    seqs_to_write += delete_extra(haplotype1_dict[gene],haplotype2_dict[gene],ref_dict[gene],phaseBlock)
        else:
            if gene in ref_dict:
                seqs_to_write += [SeqRecord(ref_dict[gene].seq,id=gene,description='')]
    SeqIO.write(seqs_to_write,args.output,'fasta')
    if args.block:
        with open(args.block,'w') as outfile:
            outfile.write("\n".join(phase_report))
            

if __name__ == "__main__":main()
            




#For sites in the longest block, replace sequences in iupac sequence with phased sequence.

#Use the intron/exon extractor as before

