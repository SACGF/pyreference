'''
Created on 19Jan.,2018

TODO: ReferenceArgs to handle changing etc




@author: dlawrence
'''

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
    


if __name__ == '__main__':
    main()