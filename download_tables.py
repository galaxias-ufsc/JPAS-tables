'''
Created on Oct 30, 2019

@author: andre
'''

from pyvo.dal import TAPService
from getpass import getpass
from astropy import log
from os import path

from cefca_tap import CEFCA_authenticate, quiet_query
from util import fix_names, convert_dtype    
    
def download_table(service, table_name, filename, overwrite=False, maxrec=100000):
    if path.exists(filename) and not overwrite:
        log.debug('File %s already exists, skipping download.')
        return
    result = quiet_query(service, 'select * from minijpas.%s' % table_name, maxrec=maxrec)
    t = result.to_table()
    fix_names(t)
    convert_dtype(t)
    del t.meta['description']
    t.write(filename, overwrite=overwrite)
    return t

log.setLevel('INFO')
tables_dir = 'tables'
maxrec = 10

table_list = ['filter',
              'magabdualobj',
              'photozlephare',
              'magabsingleobj',
              'tileimage',
              'xmatch_jplus_dr1',
              'xmatch_sdss_dr12',
              'xmatch_deep2_spec',
              'xmatch_alhambra',
              ]


login = input('Login: ')
password = getpass('Password: ')
auth = CEFCA_authenticate(login, password)
service = TAPService('https://archive.cefca.es/catalogues/vo/tap/minijpas-idr201910', auth)

for tab in table_list:
    filename = path.join(tables_dir, '%s.fits' % tab)
    log.info('Downloading table %s to %s.' % (tab, filename))
    t = download_table(service, tab, filename, overwrite=True, maxrec=10)

