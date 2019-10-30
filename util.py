'''
Created on Oct 30, 2019

@author: andre
'''

from astropy import log
import numpy as np


def fix_names(t):
    '''
    Make the table columns names simpler, and transform length=1 columns to scalars.
    '''
    for name in t.colnames:
        field_name = name.split('.')[1].lower()
        t[name].name = field_name

def convert_dtype(t):
    '''
    TAP result tables have columns with weird dtype.
    Transform the table to a more suitable format for FITS tables.
    '''
    obj_dtype = np.dtype('O')
    for d in t.dtype.descr:
        if len(d) == 3:
            name, dtype, _ = d
            t[name] = t[name].squeeze()
        elif len(d) == 2:
            name, dtype = d
        else:
            raise Exception('Unknown column description: %s' % d)
        if dtype == obj_dtype:
            log.debug('Found object type column: %s' % name)
            data = obj_col_to_array(t[name])
            t[name] = data


def obj_col_to_array(c):
    '''
    Transform a colum of dtype object, each row containing
    an array, to a 2-d array.
    '''
    # Detect the object size for initialization.
    max_len_x = 0
    for x in c:
        if len(x) > max_len_x:
            max_len_x = len(x)
            break

    if isinstance(x, bytes):
        dtype = np.dtype('|S20')
        shape = (len(c),)
        data = np.zeros(shape, dtype=dtype)
        for i, x in enumerate(c):
            if len(x) == 0:
                continue
            data[i] = x
    elif isinstance(x, np.ndarray):
        if max_len_x == 0:
            raise Exception('Column %s has all elements empty.' % c.name)
        dtype = x.dtype
        shape = (len(c), max_len_x)
        data = np.zeros(shape, dtype=dtype)
        for i, x in enumerate(c):
            data[i, :len(x)] = x
    return data
