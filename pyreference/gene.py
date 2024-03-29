"""
Created on 19Jan.,2018

@author: dlawrence
"""

from lazy import lazy
from pyreference.genomic_region import GenomicRegion
from pyreference.transcript import Transcript
import sys

try:
    _big_int = sys.maxsize  # Python 3
except AttributeError:
    _big_int = sys.maxint  # Python 2

def min_transcript_key(t):
    # We want the MAX length - and MIN ID, so sort by min but use maxint-length
    # We also want NM_007041 (len 2209) over NM_001001976 (len 2209)
    # Which is annoyingly zero padded - so use smallest ID length, then only if equal do alpha sort
    return _big_int - t.length, len(t.get_id()), t.get_id()


def min_canonical_tag(t):
    # we use 'not in' as False < True (so will get minimum)
    CANONICAL_TAGS = ["MANE Select", "MANE_Select", "RefSeq Select", "Ensembl Select"]
    return tuple([x not in t.tags for x in CANONICAL_TAGS])


class Gene(GenomicRegion):
    """ Gene (which could contain multiple transcripts) """

    @property
    def name(self):
        return self.get_gene_name()

    @property
    def description(self):
        return self._dict.get("description")

    @property
    def biotype(self):
        return self._dict.get("biotype")

    @property
    def summary(self):
        return self._dict.get("summary")

    @property
    def map_location(self):
        return self._dict.get("map_location")

    def get_gene_name(self):
        return self._dict["gene_symbol"]

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
        """ Returns longest coding transcript if gene is coding, otherwise longest transcript
            Sort transcript ID alphabetically if equal length """

        methods = {
            "tags": self.get_canonical_transcript_from_tags,
            "longest_coding": self.get_longest_coding_transcript,
            "longest": self.get_longest_transcript,
        }

        transcript = None
        for rt_method in self.reference.representative_transcript_list:
            func = methods[rt_method]
            transcript = func()
            if transcript:
                return transcript
        return transcript

    def get_canonical_transcript_from_tags(self):
        """ Using the GTF tag (eg 'MANE_select') """
        transcripts = self.transcripts
        transcripts = filter(lambda t: t.tags, transcripts)
        canonical_transcript = None
        if transcripts:
            canonical_transcript = min(transcripts, key=min_canonical_tag)
        return canonical_transcript

    def get_representative_transcript(self):
        return self.representative_transcript

    def get_longest_coding_transcript(self):
        return self.get_longest_transcript(coding_only=True)

    def get_longest_transcript(self, coding_only=False):
        transcripts = self.transcripts
        if coding_only:
            transcripts = filter(lambda t: t.is_coding, transcripts)
        
        longest_transcript = None
        if transcripts:
            longest_transcript = min(transcripts, key=min_transcript_key)
        return longest_transcript

    def __repr__(self):
        return "%s (%s) %d transcripts" % (self.get_gene_name(), self.accession_id, len(self.transcripts))
        