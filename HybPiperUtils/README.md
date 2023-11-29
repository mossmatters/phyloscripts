# HybPiper Utils

Assorted shell and Python scripts to work with data from HybPiper output.

### Software

HybPiper: [www.github.com/mossmatters/HybPiper](www.github.com/mossmatters/HybPiper)

Most of the scripts require Python 3.0 or higher and `biopython` which can be installed using `conda`.


#### `fasta_merge.py`

Script to merge (concatenate) alignments from different alignments for phylogenetic analysis.

```
usage: fasta_merge.py [-h] [--fastafiles FASTAFILES [FASTAFILES ...]]
                      [--filelist FILELIST] [--raxml {DNA,WAG,JTT,CODON}]

This script will take a list of FASTA files and concatenate them for use in
phylogenetic inference. The sequence headers (up until the first space) must be identical
in each individual FASTA file.

Individual gene sequences should be aligned prior to running this script!

This script requires BioPython to read/write FASTA sequences.

optional arguments:
  -h, --help            show this help message and exit
  --fastafiles FASTAFILES [FASTAFILES ...]
                        List of Fasta Files. Can use wildcard on Linux/Mac systems
  --filelist FILELIST   File containing list of Fasta files. Alternative to --fastalist
  --raxml {DNA,WAG,JTT,CODON}
                        Create a partition file 'partitions.raxml' intended for raxml in the current directory. For amino acid sequences, select the substitution model. To specify a separate model for 1st/2nd vs. 3rd codon positions, select CODON.
```

### `filter_by_length.py`

A script for filtering gene-sample combinations based on length filters for each gene. Two filters are available: minimum length, and percentage of mean length of the targets. The script will assume you have already run `hybpiper stats` and `hybpiper retrieve_sequences` as input.

As of HybPiper 2.1.6, the `hybpiper retrieve_sequences` command could only filter sequences at a "whole project" level - for example, removing a sample if fewer than 20% of genes were recovered.

Suggested Workflow:

1. Run `hybpiper stats` to generate the `stats.tsv` and `lengths.tsv` files
2. Run `hybpiper retrieve_sequences` to create a folder of FASTA sequences
3. Run this script to create new FASTA files based on the per-gene filters. 
    Also writes to standard output the `denylist` by gene, redirect this to save to a file.
    
The FASTA sequences will expect to have the naming scheme of HybPiper:

-    geneName.FNA for nucleotide exon files
-    geneName.FAA for amino acid files
    
The geneNames will be taken from either the `hybpiper stats` file (`--lengthfile`) or a supplied list of gene sample combinations (`--denylist`, also produced by running this script). There are two filters, `--length_filter` for the minimum length to accept a sequence (for all genes) and `--percent_filter` for a fraction of the mean length determined from the `seq_lengths.tsv` file for each gene. For example:

```
python filter_by_length.py --lengthfile ../seq_lengths.tsv --extension FNA --percent_filter 0.1 > denylist.txt  
```
    
If you wish to filter intron or supercontig sequences, run a second time with the --denylist flag to skip the filtering based on lengths.
