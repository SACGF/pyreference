from lazy import lazy

from pyreference.genomic_region import GenomicRegion


class Transcript(GenomicRegion):
    def __init__(self, *args, **kwargs):
        gene = kwargs.pop("gene", None)
        super(Transcript, self).__init__(*args, **kwargs)
        self.gene = gene
        

    @lazy
    def is_coding(self):
        feature_types = set(self._dict["exons_by_type"])
        return {"CDS", "start_codon", "stop_codon"} & feature_types


    def get_representative_transcript(self):
        return self


    def get_features_length(self, feature_type):
        length = 0
        for feature in self.get_features_in_stranded_order(feature_type):
            length += feature["end"] - feature["start"] 
        return length


    def get_features_in_stranded_order(self, feature_type):
        '''features returned sorted 5' -> 3' '''

        is_reversed = self._dict["strand"] == '-'
        if is_reversed:
            stranded_start = "end"
        else:
            stranded_start = "start"
            
        features_by_type = self._dict["exons_by_type"]
        
        return sorted(features_by_type, key=lambda x : x[stranded_start], reverse=is_reversed)

    @property
    def length(self):
        return self.get_features_length("exon")
    
    
    def get_sequence_from_features(self, feature_type):
        features = self.get_features_in_stranded_order(feature_type)
        return self.reference.get_sequence_from_features(features)

    
    
    def get_coding_sequence(self):
        ''' Warning: There are frame shift issues not handled here.
            Do not naively turn this into a protein - better to use existing databases '''
        return self.get_sequence_from_features("CDS")

    def get_5putr_sequence(self):
        return self.get_sequence_from_features("5PUTR")

    def get_3putr_sequence(self):
        return self.get_sequence_from_features("3PUTR")

    