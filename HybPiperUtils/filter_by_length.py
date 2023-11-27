import os,sys,argparse

from Bio import SeqIO

helptext ='''This script will filter output from HybPiper based on the output of hybpiper retrieve_sequences

As of HybPiper version 2.1.6, hybpiper retrieve_sequences only supports filtering based
on project-wide thresholds (i.e. number of total genes recovered). This script will allow 
filtering based on individual genes and the 

1. Run hybpiper stats to generate the stats.tsv and lengths.tsv files
2. Run hybpiper retrieve_sequences to create a folder of FASTA sequences
3. Run this script to create new FASTA files based on the per-gene filters. 
    Also makes a text file summarizing the denylist by gene.
    
The FASTA sequences will expect to have the naming scheme of HybPiper:
    geneName.FNA for nucleotide exon files
    geneName.FAA for amino acid files
    
The geneNames will be taken from either the hybpiper stats file (--lengthfile) or a supplied
    list of gene sample combinations (--denylist, also produced by running this script)
    
If you wish to filter intron or supercontig sequences, run again with the --denylist flag
    to skip the filtering based on lengths.
'''



def parse_seqlens(seqlens_fn):
    '''Takes the file name for the seqlengths output of hybpiper stats and returns:
    - a list of sample names
    - a dictionary for each gene containing:
        * the name of the gene as the dict key
        * "mean length":integer
        * "sample_lengths":[a list of sample lengths in same order as sample names list]'''
    
    sample_names = []
    gene_lengths = {}
    
    seqlens = open(seqlens_fn)
    genenames = seqlens.readline().rstrip().split("\t")[1:]
    meanlens = seqlens.readline().rstrip().split("\t")[1:]
    for geneNum in range(len(genenames)):
        gene_lengths[genenames[geneNum]] = {"mean_length":float(meanlens[geneNum]),"sample_lengths":[]}
    for line in seqlens:
        line = line.rstrip().split("\t")
        sampleName = line.pop(0)
        sample_names.append(sampleName)
        for geneNum in range(len(genenames)):
            gene_lengths[genenames[geneNum]]["sample_lengths"].append(float(line[geneNum]))
    
    return sample_names,gene_lengths


def parse_denylist(denylist_fn):
    '''parses the text file at denylist_fn and returns a dict with the geneName:[samplelist] pairs'''
    deny_dict = {}
    total_deny = 0
    for line in open(denylist_fn):
        line = line.rstrip().split("\t")
        samples = line[1].split(",")
        total_deny += len(samples)
        deny_dict[line[0]] = samples
    sys.stderr.write(f"Found {len(samples)} total samples at {len(deny_dict)} genes in the denylist {denylist_fn}")
    return deny_dict


def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--denylist",help="Text file containing gene-sample combinations to omit. \n The format of the file should be one gene per line, a tab, \n and then a comma-delimited list of samples to disallow: \n    gene[tab]sample,sample,sample ",default=None)
    parser.add_argument("--lengthfile",help="Output of hybpiper stats, with list of genes in first row, \n mean target lengths in second row, and sample recovery in other rows.")
    parser.add_argument("--extension",help="File extension for all FASTA files to filter in current directory. \n For example, the amino acid output of HybPiper would be: FAA",default="FNA")
    parser.add_argument("--length_filter",help="Minimum length to allow a sequence \n in nucleotides for DNA or amino acids for protein sequences",default=0,type=int)
    parser.add_argument("--percent_filter",help="Minimum fraction (between 0 and 1) of the mean target length to allow a sequence for a gene. \n Lengths taken from HybPiper stats file.",default=0,type=float)
    
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    if args.denylist:
        deny_dict = parse_denylist(denylist_fn)
    else:
        sample_names,gene_lengths = parse_seqlens(args.lengthfile)
        print(gene_lengths)
    
    

if __name__ == "__main__":main()