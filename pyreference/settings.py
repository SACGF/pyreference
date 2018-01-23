'''
Created on 23Jan.,2018

@author: dlawrence
'''


# Change this when you introduce breaking changes
PYREFERENCE_JSON_VERSION = 3
PYREFERENCE_JSON_VERSION_KEY = "pyreference_json_version"

CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}

# We want our features to match PyFasta so they can be used without conversion
CHROM = "chr"
START = "start"
END = "stop"
STRAND = "strand"

# Other 
IS_CODING = "is_coding"


BEST_REGION_TYPE_ORDER = ["coding", "5PUTR", "3PUTR", "non coding", "intron"]

