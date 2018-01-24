'''
Created on 22Jan.,2018

TODO: Link to wiki page about config file format

@author: dlawrence
'''

import os
from six import reraise, raise_from
from six.moves import configparser #@UnresolvedImport
import sys


def load_params_from_config(build=None, config=None):
    DEFAULT_PYREFERENCE_CFG = "~/pyreference.cfg"
    
    if config is None:
        config = os.path.expanduser(DEFAULT_PYREFERENCE_CFG)
        if not os.path.exists(config):
            msg = "No default config file: '%s' found. Please create one or pass another with the 'pyreference_cfg' parameter." % config
            raise OSError(msg)
    else:
        if not os.path.exists(config):
            msg = "Passed config file '%s' not found." % config
            raise OSError(msg)

    params = {}

    try:
        defaults = {'genes_json' : None,
                    'trna_json' : None,
                    'mature_mir_sequence_fasta' : None,
                    'genome_sequence_fasta' : None,}
        cfg = configparser.SafeConfigParser(defaults=defaults)
        cfg.read(config)

        if build is None:
            try:
                build = cfg.get("global", "default_build")
            except configparser.NoOptionError as noe:
                msg = "No build parameter passed, and no default_build set in [global] section."
                raise_from(configparser.NoOptionError(msg), noe)
                
        for k in defaults.keys():
            params[k] = cfg.get(build, k)
    except:
        msg ="Problem parsing config file '%s':" % config
        traceback = sys.exc_info()[2]
        reraise(configparser.NoOptionError, msg, traceback)
    
    return params