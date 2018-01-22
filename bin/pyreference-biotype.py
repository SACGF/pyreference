'''
Created on 19Jan.,2018

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



if __name__ == '__main__':
    main()