"""
Created on 22Jan.,2018

TODO: Link to wiki page about config file format

@author: dlawrence
"""

import os
from six import raise_from
from six.moves import configparser  # @UnresolvedImport

try:
    from ConfigParser import SafeConfigParser as ConfigParser
except ImportError:
    from configparser import ConfigParser


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

    GLOBAL_FLAGS = ["use_gzip_open", "stranded"]
    params = {}

    defaults = {
        'genome_accession': None,
        'genes_json': None,
        'trna_json': None,
        'mature_mir_sequence_fasta': None,
        'genome_sequence_fasta': None,
        "genome_sequence_lookup": None,
    }
    cfg = ConfigParser(allow_no_value=True, defaults=defaults)
    cfg.read(config)

    if build is None:
        try:
            build = cfg.get("global", "default_build")
        except configparser.NoOptionError as noe:
            msg = "No build parameter passed, and no default_build set in [global] section of config file '%s'" % config
            raise_from(configparser.NoOptionError(msg), noe)

    if not cfg.has_section(build):
        msg_params = {"build": build, "config": config}
        msg = "Build='%(build)s', no section [%(build)s] in config file '%(config)s" % msg_params
        raise ValueError(msg)

    for f in GLOBAL_FLAGS:
        try:
            params[f] = cfg.getboolean("global", f)
        except configparser.NoOptionError as noe:
            pass

    params["build"] = build
    params["config"] = config

    for k in defaults.keys():
        params[k] = cfg.get(build, k)

    return params
