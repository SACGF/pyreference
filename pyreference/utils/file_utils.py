'''
Created on 22Jan.,2018

@author: dlawrence
'''
import os


try:
    from pathlib import Path
except (ImportError,AttributeError):
    from pathlib2 import Path

def name_from_file_name(file_name):
    return Path(file_name).name

def stem_from_file_name(file_name):
    return Path(file_name).stem

def mk_path(path):
    if path and not os.path.exists(path):
        os.makedirs(path)

def mk_path_for_file(f):
    mk_path(os.path.dirname(f))
    
def file_or_file_name(f, mode='r'):
    if isinstance(f, basestring): # Works on unicode
        if 'w' in mode: # Create path if writing
            mk_path_for_file(f)
            
        return open(f, mode)
    elif isinstance(f, file):
        return f # Already a File object
    else:
        raise ValueError("'%s' (%s) not a file or string" % (f, type(f)))
