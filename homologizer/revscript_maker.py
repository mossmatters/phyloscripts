

# This script will generate the RevBayes scripts necessary to run Homologizer on sets of genes.
# 
# Matt Johnson, Texas Tech University


######### TO DO / CONSIDERATIONS
#
# RevBayes can read the label-swap from a txt file

# What will RevBayes do when a taxon that has homeologs is missing?
# There's code to fill in missing taxa but I suppose it will fill in both alleles
# This SHOULD result in a 50/50 probability for that gene
# Will be accommodated by the script that does the label switching, need to remember that
#    maybe the taxon isn't present
# I'm also not sure how many genes will only have one homeolog (and therefore no _h1 or _h2)
# This way, should be able to use one taxon-swap file for all genes
# First version: just use genes with complete sampling.


# Also need to get it to read the correct alignments
# Is it possible to pass command line arguments to revbayes scripts?
#    Guessing not, so Python can write a header for the script call the appropriate
#    
# When doing subsets, how do we know that the phase is the same from subset to subset?
#   Could do one final analysis with one gene picked from each of the subsets

import argparse,pathlib,os
from Bio import AlignIO,SeqIO
from Bio.Alphabet import generic_dna,generic_protein
from Bio.Nexus import Nexus
from io import StringIO


helptext = '''
###### REQUIREMENTS
#
#   Python > 3.5
#   Biopython

###### INPUTS
#
#   Aligned sequence files in a format readable by BioPython
#   Max number of genes to include per RevBayes run

##### OUTPUTS
#   For each run of RevBayes:
#       a) a text file containing the locations of the gene alignments in the subset
#       b) a RevBayes script to run a subset of genes
#   A text file containing the names of the RevBayes scripts (could be used for SGE array jobs)
'''

def convert_to_nexus(alignment_fn,genedir,input_type="fasta",seqtype = 'dna'):
    if seqtype == 'dna':
        seqs_to_write = [seq for seq in SeqIO.parse(os.path.join(genedir,alignment_fn),input_type, alphabet=generic_dna)]
    else:
        seqs_to_write = [seq for seq in SeqIO.parse(os.path.join(genedir,alignment_fn),input_type, alphabet=generic_protein)]
    new_fn = ".".join(alignment_fn.split(".")[:-1]) + ".nexus"    
    output = StringIO()
    SeqIO.write(seqs_to_write, output, 'nexus')
    p = Nexus.Nexus()
    p.read(output.getvalue())
    p.write_nexus_data(os.path.join("nexusfiles",new_fn), interleave=False)



parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("genedir",help="directory containing all of the gene alignments. No other files should be in this directory.")
parser.add_argument("swaplist",help="file containing labels to swap, two per sample, separated by a comma")
parser.add_argument("--numgenes","-n",help="number of genes to include in each RevBayes job",type=int,default=25)
parser.add_argument("--alignfiletype","-a",help="Alignment file type. Must be one used by BioPython",default="fasta")
parser.add_argument("revbayestemplate",help="Template for RevBayes. This script will prepend with alignment and swap info")
args = parser.parse_args()



alignments = os.listdir(args.genedir)

######## CONVERT TO NEXUS

pathlib.Path("nexusfiles").mkdir(parents=True, exist_ok=True)

if args.alignfiletype != "nexus":
    for fn in alignments:
        convert_to_nexus(fn,args.genedir)
    
split_gene_list = [alignments[x:x+args.numgenes] for x in range(0, len(alignments), args.numgenes)]
print("Will generate {} RevBayes scripts with a max of {} genes".format(len(split_gene_list),args.numgenes))

######## WRITE REVBAYES SCRIPTS
revbayes_header = '''
####ADD THESE FROM PYTHON SCRIPT
output_file = "alleles.{}"
label_swap = readTable("{}",delimiter=",")
alignments = readTable("{}")
#######
'''

pathlib.Path("revbayes_scripts").mkdir(parents=True, exist_ok=True)
revbayes_text = open(args.revbayestemplate).read()

for gs in range(len(split_gene_list)):
    genelist_file = "revbayes_scripts/genelist.{}.txt".format(gs)
    with open(genelist_file,'w') as outfile:
        for gene in split_gene_list[gs]:
            outfile.write(args.genedir + "/" + gene+"\n")
            
    with open("revbayes_scripts/homeolog_phase.{}.Rev".format(gs),'w') as revbayes_out:
        revbayes_out.write(revbayes_header.format(gs,args.swaplist,genelist_file))
        revbayes_out.write(revbayes_text)





