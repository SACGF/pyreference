"""

ACRF Cancer Genomics Facility

Created on 23Jan.,2018

@see http://www.mirbase.org/

@author: dlawrence
"""

from Bio.Seq import Seq

def rna_to_dna(rna_sequence):
    rna = Seq(rna_sequence)
    return rna.back_transcribe()

def rna_to_target(rna_sequence):
    dna = rna_to_dna(rna_sequence)
    return str(dna.reverse_complement())

class MiRNA(object):
    def __init__(self, name, ref_sequences):
        self.name = name
        self.mature = ref_sequences[self.name]

    def get_8mer_target(self):
        return rna_to_target(self.mature[:8])

    # Nucleotides 4-15
    def get_central_paired_target(self):
        return rna_to_target(self.mature[3:15])

    def get_mature_dna(self):
        return str(rna_to_dna(self.mature))

