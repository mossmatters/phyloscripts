# Homologizer

Uses the setHomeologPhase function in [RevBayes](revbayes.github.io) to adjust the labels in datasets that contain homeologs. Once adjusted, the data from many genes can be analyzed together, for example using a concatenated supermatrix or summary species tree analysis. 

### Step 0: Input files

Gene alignments that include samples that have two homeolog sequences. 
The sequences within each gene should be phased. 
For target capture data, see the `alleles_workflow` and `haplonerate` methods in this same [Phyloscripts](https://github.com/mossmatters/phyloscripts) repository for tips on generating phased haplotypes within a gene sequence. The gene alignments can be in FASTA or NEXUS format, but must end with a regular suffix. For example, all fasta files must end with `.fasta` or `.fa`. Place all of the alignment files in a directory. **Do not place any other files in that directory**.

Users must also prepare a text file containing the naming scheme for sequences. For example:

	Entosthodon-hungaricus-3838_h1	Entosthodon-hungaricus-3838-h2
	Physcomitrium-immersum-3176_h1	Physcomitrium-immersum-3176_h2
	Physcomitrium-pyriforme-3410_h1	Physcomitrium-pyriforme-3410_h2

These labels will be used during the RevBayes script to swap labels for the polyploid individuals.

**In this version, only two sequences per individual are supported**

### Step 1: Making RevBayes scripts

Given a set of gene alignments, we will need to prepare RevBayes files that have the appropriate alleles to switch. The paired labels are taken from the label swap file, and is added to a basic template for RevBayes (taken from Will Freyman's version).

**For now, only genes with no missing labels are accepted**

For ease of analysis, all the genes are split up into chunks of genes (by default, max 25 genes). This means that for 250 loci, there will be 10 RevBayes scripts generated.

Sequences are converted NEXUS for RevBayes.

### Step 2: Run RevBayes on sets of genes
   
Issues:

- MPIrun version of RevBayes crashes, probably a memory error
- Single-thread version of RevBayes grabs all RAM on the machine (256 GB on one machine!)

### Step 3: Summarize RevBayes output


The script `swap_labels.py` reads the RevBayes log files and calculates the posterior probability of label swapping. A threshold can be picked (set to 95% by default), if the swapping PP is below this, both sequences are deleted from that sample. The output 

# TODO:

- Incorporate a burnin to the `swap_labels.py` script
- Get RevBayes running more efficiently
- Accomodate missing labels
- How to summarize across different chunks of genes. Each chunk will be homologized, but no guarantee of this across chunks.