'''
Created on 19Jan.,2018

@author: dlawrence
'''

from lazy import lazy

from pyreference.genomic_region import GenomicRegion
from pyreference.transcript import Transcript


class Gene(GenomicRegion):
    ''' Gene (which could contain multiple transcripts) '''

    @property
    def name(self):
        return self.get_gene_name()

    def get_gene_name(self):
        return self._dict["name"]

    @property
    def biotype(self):
        return '/'.join(sorted(self.get_biotypes()))

    def get_biotypes(self):
        return self._dict["biotype"]

    @lazy
    def transcripts(self):
        transcripts = []
        
        for transcript_id in self._dict["transcripts"]:
            td = self.reference.get_transcript_dict(transcript_id)  
            transcript = Transcript(self.reference, transcript_id, td, gene=self)
            transcripts.append(transcript)
        return transcripts


    @lazy
    def is_coding(self):
        return any(t.is_coding for t in self.transcripts)


    @lazy
    def representative_transcript(self):
        ''' Returns longest coding transcript if gene is coding, otherwise longest transcript '''
        
        transcript = self.get_longest_coding_transcript()
        if transcript == None:
            transcript = self.get_longest_transcript()
        return transcript
    

    def get_representative_transcript(self):
        ''' Representative transcript = longest coding transcript if gene is coding, otherwise longest transcript
            Returns one Transcript instance or None
        '''
        return self.representative_transcript
        

    def get_longest_coding_transcript(self):
        return self.get_longest_transcript(coding_only=True)

    
    def get_longest_transcript(self, coding_only=False):
        transcripts = self.transcripts
        if coding_only:
            transcripts = filter(transcripts, lambda t : t.is_coding)
        return sorted(transcripts, key=lambda t : t.length, reversed=True)


    def __repr__(self):
        return "%s (%s) %d transcripts" % (self.get_gene_name(), self.accession_id, len(self.transcripts))
        