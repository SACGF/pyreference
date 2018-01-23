'''
Created on 19Jan.,2018

TODO: ReferenceArgs to handle changing etc




@author: dlawrence
'''

from inspect import getsourcefile
import os
from os.path import abspath

from pyreference import Reference


def main():
    # Use this as a test platform to load reference
    reference = Reference()

    msn = reference.get_gene_by_name("MSN")
    print(msn) 
    print("promoter: %s" % msn.get_promoter_sequence())


    gata2 = reference.get_gene_by_name("GATA2")
    print(gata2) 
    print("promoter: %s" % gata2.get_promoter_sequence())


    print("----------")
    
    for (gene_id, gene) in reference.protein_coding_genes.items():
        print("gene_id=%s, gene=%s" % (gene_id, gene)) 
    

    this_file_dir = os.path.dirname(abspath(getsourcefile(lambda:0)))
    reference_dir = os.path.join(this_file_dir, "..", "pyreference", "tests", "reference")
    
    genes_json = os.path.join(reference_dir, "hg19_chrY_300kb_genes.gtf.json.gz")
    genome_sequence_fasta = os.path.join(reference_dir, "hg19_chrY_300kb.fa")
    reference = Reference(genes_json=genes_json,
                          genome_sequence_fasta=genome_sequence_fasta)

    t = reference.get_transcript_by_id("NR_028057_2")
    print("%s is_coding: %s" % (t, t.is_coding))





if __name__ == '__main__':
    main()