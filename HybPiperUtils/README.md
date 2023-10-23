# HybPiper Utils

Assorted shell and Python scripts to work with data from HybPiper output.

### Software

HybPiper: www.github.com/mossmatters/HybPiper

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