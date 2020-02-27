#Script to swap labels in alignments based on the results of RevBayes


helptext = '''
###### Input
#
# List of genes (from first script)
# RevBayes output log files

###### Options
#
# Swap threshold- based on posterior distribution, below threshold sample is deleted (default 95%)
# Burnin percentage (default 10%)
# Output alignment file type (default FASTA)

##### Output
#
# Directory containing alignments of specified type with the labels swapped.

'''

import argparse,pathlib,os
import pandas as pd
from Bio import SeqIO, AlignIO

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('genelist',help = "File containing list of gene alignments.")
parser.add_argument('initalAlignmentDir',help="Directory containing RevBayes logs for each alignment")
parser.add_argument('--swappct','-s',help="Posterior probability to do swapping. Below this, both alleles from the sample are deleted.",default=95)
parser.add_argument('--burnin','-b',help="Percentage of RevBayes to discard as burnin",default=10)
parser.add_argument('--outfiletype','-f',help="Alignment file type to output",default='fasta')
parser.add_argument('--infiletype','-i',help="Alignment file type to output",default='fasta')

args = parser.parse_args()


###### READ IN GENELIST

genelist = [x.rstrip() for x in open(args.genelist)]
subset_num = args.genelist.split(".")[1]
outputdir = args.initalAlignmentDir+"_swapped"
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)

for gene in genelist:
###### READ IN REVBAYES LOGS
    logfilepath = "alleles.{}_{}_phase.log".format(subset_num,gene)
    gene_logfile = pd.read_csv(logfilepath,header=0,index_col=0,sep="\t") 
    dim_logfile = gene_logfile.shape
    burnin = int(dim_logfile[0]/args.burnin)
    gene_swap_dict = {}
    for h in range(dim_logfile[1]):
        pp = gene_logfile.iloc[burnin:,h].value_counts()/dim_logfile[0] * 100
        if pp[0] <  args.swappct:
            #None will indicate the sequence should be skipped when re-writing
            gene_swap_dict[pp.name] = None
            print("PP for {} in gene {} was {}".format(pp.name,gene,pp[0]))
        else:
            gene_swap_dict[pp.name] = pp.index[0]
    with open(os.path.join(outputdir,os.path.split(gene)[1]),'w') as outfile:
        for seq in SeqIO.parse(gene,args.infiletype):
            if seq.id in gene_swap_dict:
                if gene_swap_dict[seq.id]:
                    seq.id = gene_swap_dict[seq.id]
                    SeqIO.write(seq,outfile,args.outfiletype)
            else:
                SeqIO.write(seq,outfile,args.outfiletype)

##### READ IN ALIGNMENT AND SWAP LABELS

