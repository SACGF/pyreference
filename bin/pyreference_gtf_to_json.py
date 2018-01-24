#!/usr/bin/env python

from __future__ import print_function, absolute_import

import HTSeq
from collections import defaultdict
import gzip
import json
import os
from pyreference.settings import CHROM, START, END, STRAND, IS_CODING, CODING_FEATURES, \
    PYREFERENCE_JSON_VERSION_KEY, PYREFERENCE_JSON_VERSION
from pyreference.utils.file_utils import name_from_file_name
import sys


class SetEncoder(json.JSONEncoder):
    ''' Dump set as list, from: https://stackoverflow.com/a/8230505/295724 '''
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)

def get_biotype_from_transcript_id(transcript_id):
    biotypes_by_transcript_id_start = {"NM_" : "protein_coding", "NR_" : "non_coding"}
    for (start, biotype) in biotypes_by_transcript_id_start.items():
        if transcript_id.startswith(start):
            return biotype

    if "tRNA" in transcript_id:
        return "tRNA"
    return "N/A"

def update_extents(genomic_region_dict, feature):
    start = genomic_region_dict[START]
    if feature.iv.start < start:
        genomic_region_dict[START] = feature.iv.start  
    
    end = genomic_region_dict[END]
    if feature.iv.end > end:
        genomic_region_dict[END] = feature.iv.end


def add_UTR_features(transcripts_by_id, transcript_cds_by_id):
    ''' Add 5PUTR/3PUTR features to coding transcripts '''

    for transcript_id, transcript in transcripts_by_id.items():
        cds_extent = transcript_cds_by_id.get(transcript_id)
        if cds_extent:
            transcript[IS_CODING] = 1
            features_by_type = transcript["features_by_type"]
            
            (left, right) = ("5PUTR", "3PUTR")
            if transcript[STRAND] == '-': # Switch
                (left, right) = (right, left)
            
            cds_start = cds_extent[START]
            cds_end = cds_extent[END]
            
            for exon in features_by_type["exon"]:
                exon_start = exon[START] 
                exon_end = exon[END] 
                
                if exon_start < cds_start:
                    end_non_coding = min(cds_start, exon_end)
                    utr_feature = {START : exon_start,
                                   END : end_non_coding}
                    features_by_type[left].append(utr_feature)
        
                if exon_end > cds_end:
                    start_non_coding = max(cds_end, exon_start)
                    utr_feature = {START : start_non_coding,
                                   END : exon_end} 
                    features_by_type[right].append(utr_feature)


def main():
    if len(sys.argv) != 2:
        script_name = os.path.basename(__file__)
        print("Usage: %s reference.gtf" % script_name, file=sys.stderr)
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    
    genes_by_id = {}
    transcripts_by_id = {}
    # Store CDS in separate dict as we don't need to write as JSON
    transcript_cds_by_id = {} 
    gene_id_by_name = {}
    gene_ids_by_biotype = defaultdict(set)
    
    
    for feature in HTSeq.GFF_Reader(gtf_file):
        gene_id = feature.attr["gene_id"]
        
        gene_name = feature.attr["gene_name"]
        gene_id_by_name[gene_name] = gene_id # TODO: Check for dupes?

        gene = genes_by_id.get(gene_id)
        if gene is None:
            gene = {"name" : gene_name,
                    "transcripts" : set(),
                    "biotype" : set(),
                    CHROM : feature.iv.chrom,
                    START : feature.iv.start,
                    END : feature.iv.end,
                    STRAND : feature.iv.strand,}
        
            genes_by_id[gene_id] = gene
        else:
            update_extents(gene, feature)
    
        transcript_id = feature.attr["transcript_id"]
        gene["transcripts"].add(transcript_id)
        transcript = transcripts_by_id.get(transcript_id)
        if transcript is None:
            transcript = {"features_by_type" : defaultdict(list),
                          "biotype" : set(),
                          CHROM : feature.iv.chrom,
                          START : feature.iv.start,
                          END : feature.iv.end,
                          STRAND : feature.iv.strand,
                          IS_CODING : 0,}
            transcripts_by_id[transcript_id] = transcript
        else:
            update_extents(transcript, feature)

        # No need to store chrom/strand for each feature, will use transcript         
        feature_dict = {START : feature.iv.start,
                        END : feature.iv.end,}
    
        transcript["features_by_type"][feature.type].append(feature_dict)
        if feature.type in CODING_FEATURES:
            cds_extent = transcript_cds_by_id.get(transcript_id)
            if cds_extent is None:
                cds_extent = {START : feature.iv.start,
                              END : feature.iv.end}
                transcript_cds_by_id[transcript_id] = cds_extent
            else:
                cds_extent[START] = min(cds_extent[START], feature.iv.start)
                cds_extent[END] = max(cds_extent[END], feature.iv.end)
        
        biotype = feature.attr.get("gene_biotype")
        if biotype is None:
            biotype = get_biotype_from_transcript_id(transcript_id)
        
        if biotype:
            gene["biotype"].add(biotype)
            transcript["biotype"].add(biotype)
        
        gene_ids_by_biotype[biotype].add(gene_id)
    
    add_UTR_features(transcripts_by_id, transcript_cds_by_id)
    
    data = {PYREFERENCE_JSON_VERSION_KEY : PYREFERENCE_JSON_VERSION,
            "genes_by_id" : genes_by_id,
            "transcripts_by_id" : transcripts_by_id,
            "gene_id_by_name" : gene_id_by_name,
            "gene_ids_by_biotype" : gene_ids_by_biotype,}
    
    genes_json_gz = name_from_file_name(gtf_file) + ".json.gz"
    with gzip.open(genes_json_gz, 'w') as outfile:
        json_str = json.dumps(data, cls=SetEncoder)
        outfile.write(json_str.encode('ascii')) 

    print("Wrote:", genes_json_gz)
    
if __name__ == '__main__':
    main()