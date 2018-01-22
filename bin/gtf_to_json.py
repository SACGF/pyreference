#!/usr/bin/env python3

from __future__ import print_function

import HTSeq
from collections import defaultdict
import gzip
import json
import sys
from pyreference.utils.file_utils import name_from_file_name

try:
    from pathlib import Path
except (ImportError,AttributeError):
    from pathlib2 import Path


class SetEncoder(json.JSONEncoder):
    ''' Dump set as list, from: https://stackoverflow.com/a/8230505/295724 '''
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)

def create_transcript():
    return {"features_by_type" : defaultdict(list),
            "biotype" : set()}


def get_biotype_from_transcript_id(transcript_id):
    biotypes_by_transcript_id_start = {"NM_" : "protein_coding", "NR_" : "non_coding"}
    for (start, biotype) in biotypes_by_transcript_id_start.items():
        if transcript_id.startswith(start):
            return biotype

    if "tRNA" in transcript_id:
        return "tRNA"
    return "N/A"


def main():
    if len(sys.argv) != 2:
        print(sys.stderr, "Usage: %s reference.gtf" % sys.argv[0])
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    
    genes_by_id = {}
    transcripts_by_id = defaultdict(create_transcript)
    gene_id_by_name = {}
    gene_ids_by_biotype = defaultdict(set)
    
    for feature in HTSeq.GFF_Reader(gtf_file):
        gene_id = feature.attr["gene_id"]
        
        gene = genes_by_id.get(gene_id)
        gene_name = feature.attr["gene_name"]
        gene_id_by_name[gene_name] = gene_id # TODO: Check for dupes?
        
        if gene is None:
            gene = {"name" : gene_name,
                    "transcripts" : set(),
                    "chrom" : feature.iv.chrom,
                    "start" : feature.iv.start,
                    "end" : feature.iv.end,
                    "strand" : feature.iv.strand,}
        
            genes_by_id[gene_id] = gene
        else:
            # Update start/end
            start = gene["start"]
            if feature.iv.start < start:
                gene["start"] = feature.iv.start  
            
            end = gene["end"]
            if feature.iv.end > end:
                gene["end"] = feature.iv.end  
    
    
        transcript_id = feature.attr["transcript_id"]
        gene["transcripts"].add(transcript_id)
        transcript = transcripts_by_id[transcript_id]
        exon_dict = {"start" : feature.iv.start,
                     "end" : feature.iv.end,
                     "num" : feature.attr["exon_number"]}
    
        transcript["features_by_type"][feature.type].append(exon_dict)
        biotype = feature.attr.get("gene_biotype")
        if biotype is None:
            biotype = get_biotype_from_transcript_id(transcript_id)
        
        if biotype:
            transcript["biotype"].add(biotype)
        
        gene_ids_by_biotype[biotype].add(gene_id)
        
    
    data = {"genes_by_id" : genes_by_id,
            "transcripts_by_id" : transcripts_by_id,
            "gene_id_by_name" : gene_id_by_name,
            "gene_ids_by_biotype" : gene_ids_by_biotype,}
    
    genes_json_gz = name_from_file_name(gtf_file) + ".json.gz"
    with gzip.open(genes_json_gz, 'w') as outfile:
        json.dump(data, outfile, cls=SetEncoder)

    print("Wrote:", genes_json_gz)
    
if __name__ == '__main__':
    main()