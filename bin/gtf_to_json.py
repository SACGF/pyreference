#!/usr/bin/env python3

from __future__ import print_function
import HTSeq
from collections import defaultdict
import json


REFERENCE_FILE = "/media/dlawrence/SDHD2/reference/ensembl_grch38_genes.gtf"

genes = {}

class SetEncoder(json.JSONEncoder):
    ''' Dump set as list, from: https://stackoverflow.com/a/8230505/295724 '''
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def main():
    def create_transcript():
        return {"exons_by_type" : defaultdict(list),
                "biotype" : set()}
    
    for feature in HTSeq.GFF_Reader(REFERENCE_FILE):
        gene_id = feature.attr["gene_id"]
        
        gene = genes.get(gene_id)
        if gene is None:
            gene = {"name" : feature.attr["gene_name"],
                    "transcripts" : defaultdict(create_transcript),
                    "start" : feature.iv.start,
                    "end" : feature.iv.end,
                    "strand" : feature.iv.strand,}
        
            genes[gene_id] = gene
        else:
            # Update start/end
            start = gene["start"]
            if feature.iv.start < start:
                gene["start"] = feature.iv.start  
            
            end = gene["end"]
            if feature.iv.end > end:
                gene["end"] = feature.iv.end  
    
    
        transcript_id = feature.attr["transcript_id"]
        transcript = gene["transcripts"][transcript_id]
        exon_dict = {"start" : feature.iv.start,
                     "end" : feature.iv.end,
                     "num" : feature.attr["exon_number"]}
    
        transcript["exons_by_type"][feature.type].append(exon_dict)
        biotype = feature.attr["gene_biotype"]
        transcript["biotype"].add(biotype)
    
    
    with open('transcripts_by_gene.json', 'w') as outfile:
        json.dump(genes, outfile, cls=SetEncoder)

    
if __name__ == '__main__':
    main()