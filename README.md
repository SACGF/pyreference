## PyReference ##

A Python library for working with reference gene annotations.

PyReference loads GTF annotations extremely rapidly, and makes it easy to write code which can be run against different genomes.

## Example ##

    import numpy as np
    import pyreference
	
	reference = pyreference.Reference()

	my_gene_ids = ["MSN", "GATA2", "ZEB1"]
	for gene in reference[my_gene_ids]:
		average_length = np.mean([t.length for t in gene.transcripts])
		print("%s average length = %.2f" % (gene, average_length))
		print(gene.iv)
		for transcript in gene.transcripts:
			if transcript.is_coding:
				threep_utr = transcript.get_3putr_sequence()
				print("%s end of 3putr: %s" % (transcript.get_id(), threep_utr[-20:]))

Outputs:

	MSN (MSN) 1 transcripts average length = 3970.00
	chrX:[64887510,64961793)/+
	NM_002444 end of 3putr: TAAAATTTAGGAAGACTTCA

	GATA2 (GATA2) 3 transcripts average length = 3367.67
	chr3:[128198264,128212030)/-
	NM_001145662 end of 3putr: AATACTTTTTGTGAATGCCC
	NM_001145661 end of 3putr: AATACTTTTTGTGAATGCCC
	NM_032638 end of 3putr: AATACTTTTTGTGAATGCCC

	ZEB1 (ZEB1) 6 transcripts average length = 6037.83
	chr10:[31608100,31818742)/+
	NM_001174093 end of 3putr: CTTCTTTTTCTATTGCCTTA
	NM_001174094 end of 3putr: CTTCTTTTTCTATTGCCTTA
	NM_030751 end of 3putr: CTTCTTTTTCTATTGCCTTA
	NM_001174096 end of 3putr: CTTCTTTTTCTATTGCCTTA
	NM_001174095 end of 3putr: CTTCTTTTTCTATTGCCTTA
	NM_001128128 end of 3putr: CTTCTTTTTCTATTGCCTTA

This takes less than 4 seconds to load via a network drive on my machine.

## Installation ##

    sudo pip install pyreference

Pre-process your GTF files to create genes.gtf.json.gz (~1/20th the size of the input GTF file)

    pyreference_gtf_to_json.py genes.gtf

Create a ~/pyreference.cfg file pointing to your references.

	[global]
	default_build=hg19

	[hg19]
	genes_json=/data/reference/hg19/genes.gtf.json.gz
	mature_mir_sequence_fasta=/data/reference/hg19/mature.fa
	genome_sequence_fasta=/data/reference/hg19/genome.fa

	[mm10]
	genes_json=/data/reference/mm10/genes.gtf.json.gz
	mature_mir_sequence_fasta=/data/reference/mm10/mature.fa
	genome_sequence_fasta=/data/reference/mm10/genome.fa


## Command line arguments ##

Substitute ArgumentParser with pyreference.ReferenceArgumentParser to add a --build option to your command line arguments. 

args.reference is now initialised to the correct build/annotation.


	from pyreference import ReferenceArgumentParser

	parser = ReferenceArgumentParser()
	parser.add("mirna_name")

	args = parser.parse_args()
	reference = args.reference.get_mirna(args.mirna_name)
	print(mir.get_8mer_target())

