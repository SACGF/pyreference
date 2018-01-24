'''
Created on 24Jan.,2018

@author: dlawrence
'''
import configargparse
from lazy import lazy
from pyreference import Reference


class ReferenceArgumentParser(configargparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ReferenceArgumentParser, self).__init__(*args, **kwargs)
        self.add("--build", env_var="BUILD", help='Use [build] section of config file.')
        self.add("--config", is_config_file=True, help='Config file used for [build] section')
        self.add("--genes_json", help='gtf.json.gz file created by PyReference gtf_to_json')
        self.add("--trna_json", help='gtf.json.gz file created by PyReference gtf_to_json')
        self.add("--mature_mir_sequence_fasta", help='MiRNA Fasta') 
        self.add("--genome_sequence_fasta", help='Genome Sequence Fasta')


    @lazy
    def reference(self):
        args = self.parse_args()
        return Reference(**args.__dict__)
            
            
    

