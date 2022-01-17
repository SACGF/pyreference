"""
Created on 22Jan.,2018

@author: dlawrence
"""
from hashlib import md5
import os
import six


try:
    from pathlib import Path  # @UnresolvedImport
except (ImportError,AttributeError):
    from pathlib2 import Path  # @UnresolvedImport


def name_from_file_name(file_name):
    """ /path/to/foo.bam => foo.bam """
    return Path(file_name).name


def stem_from_file_name(file_name, remove_gz_first=False):
    if remove_gz_first and file_name.endswith(".gz"):
        file_name = file_name[:-3]
    return Path(file_name).stem


def mk_path(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


def mk_path_for_file(f):
    mk_path(os.path.dirname(f))
    

def file_or_file_name(f, mode='r'):
    if isinstance(f, six.string_types):
        if 'w' in mode:  # Create path if writing
            mk_path_for_file(f)
            
        return open(f, mode)
    elif hasattr(f, 'read'):
        return f  # Already a File object
    else:
        raise ValueError("'%s' (%s) not a file or string" % (f, type(f)))


def file_md5sum(filename):
    m = md5()
    with open(filename, "rb") as f:
        m.update(f.read())
    return m.hexdigest()
