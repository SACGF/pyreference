'''
Created on 22Jan.,2018

@author: dlawrence
'''
import csv

from pyreference.utils.file_utils import file_or_file_name


def write_csv_dict(csv_file, headers, rows, extrasaction=None, dialect=None):
    '''
    default dialect = 'excel', other dialect option: 'excel-tab'
    These are the same optional arguments as csv.DictWriter
    headers=keys for dicts
    rows=list of dicts
    '''

    
    if extrasaction is None:
        extrasaction = "raise"
    if dialect is None:
        dialect = 'excel'
    
    f = file_or_file_name(csv_file, "wb")

    writer = csv.DictWriter(f, headers, extrasaction=extrasaction, dialect=dialect)
    writer.writerow(dict(zip(headers, headers)))
    writer.writerows(rows)
