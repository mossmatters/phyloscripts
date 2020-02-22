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
    Given a set of gene alignments, we will need to prepare RevBayes files that have the appropriate alleles to switch
    Not every sample will be present in every gene so need to account for missing labels when running RevBayes
    Need to convert FASTA to NEXUS as before.

### Step 2: Run RevBayes on sets of genes
    The release version for Linux uses a Singularity image - get this working on Quanah
    Write a QSUB script to distribute the RevBayes jobs in an array

### Step 3: Summarize RevBayes output
    The script should incorporate the bayesian posterior probability of the label IDs for each gene
        I am thinking of genes where the two E. hungaricus alleles are actually sister to each other (probably gene duplication since polyploidy event).
    Do we just kick out genes where not everything has a 100% PP for label switching?
        Or do we just cut out the taxon?
        Should summarize this so we know how many genes have “clear” label switching for each allopolyploid sample

### Step 4: Change labels in RAXML collapsed gene trees

### Step 5: Run ASTRAL using gene trees with phased labels
