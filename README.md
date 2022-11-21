## PyReference ##

[![PyPi version](https://img.shields.io/pypi/v/pyreference.svg)](https://pypi.org/project/pyreference/) [![Python versions](https://img.shields.io/pypi/pyversions/pyreference.svg)](https://pypi.org/project/pyreference/)

A Python library for working with reference gene annotations. For RefSeq/Ensembl GRCh37/GRCh38 and other species

A GTF/GFF3 can take minutes to load. We pre-process it into JSON, so it can be loaded extremely rapidly.  

PyReference makes it easy to write genomics code, which is easily run across different genomes or annotation versions.

## Example ##

    import numpy as np
    from pyreference import Reference 
    
    reference = Reference()  # uses ~/pyreference.cfg default_build

    my_gene_symbols = ["MSN", "GATA2", "ZEB1"]
    for gene in reference[my_gene_symbols]:
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

This takes 4 seconds to load on my machine.

## pyreference biotype ##

Also included is a command line tool (pyreference_biotype.py) which shows which biotypes small RNA fragments map to.

![](https://i.stack.imgur.com/Tsjr3.jpg)

## Installation ##

    sudo pip install pyreference

Then you will need to:

* [Download / Create gene annotations](https://github.com/SACGF/pyreference/wiki/genes_json_file)
* Create a [pyreference config files](https://github.com/SACGF/pyreference/wiki/pyreference_config_file)
