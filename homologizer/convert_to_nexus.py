
# Script to convert Physcomitrium/Entosthodon phased FASTA gene files into NEXUS format
# Only samples in the list will be retained
# The gene name will also be stripped when writing out
import sys
from Bio import AlignIO,SeqIO
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
from io import StringIO



samples_to_keep = ["Physcomitrium-immersum-3176_h1",
"Physcomitrium-immersum-3176_h2",
"Entosthodon-hungaricus-3838_h1",
"Entosthodon-hungaricus-3838_h2",
"Physcomitrium-pyriforme-3410_h1",
"Physcomitrium-pyriforme-3410_h2"]

alignment_fn = sys.argv[1]
geneID = alignment_fn.split(".")[0]

reduced_alignment_fn = "{}.revbayes.nexus".format(geneID)

seqs_to_write = []
with open(reduced_alignment_fn,'w') as outfile:
	for seq in SeqIO.parse(alignment_fn,'fasta', alphabet=generic_dna):
		seq.id = seq.id.replace("-{}".format(geneID),"")

		if "_h" in seq.id:
			if seq.id in samples_to_keep:
				seqs_to_write.append(seq)
		else:
			seqs_to_write.append(seq)
	#SeqIO.write(seqs_to_write,reduced_alignment_fn,'nexus')

output = StringIO()
SeqIO.write(seqs_to_write, output, 'nexus')
p = Nexus.Nexus()
p.read(output.getvalue())
p.write_nexus_data(reduced_alignment_fn, interleave=False)

#AlignIO.convert(reduced_alignment_fn,'fasta',reduced_alignment_fn.replace("fasta","nexus") ,'nexus',generic_dna,interleave=False)

