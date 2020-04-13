'''
Created on Oct 30, 2019

@author: andre
'''

import argparse
from getpass import getpass
from jpas_tap import TAPQueueManager
from astropy import log
from time import sleep
from os import path

default_service_url = 'https://archive.cefca.es/catalogues/vo/tap/minijpas-pdr201912'

def parse_args():
    parser = argparse.ArgumentParser(description='Download JPAS tables.')
    parser.add_argument('objectList', type=str, nargs='*',
                        help='Fit a list of objects, ignoring --range option.')
    parser.add_argument('--login', dest='login', help='Authentication username. Will prompt if not set.')
    parser.add_argument('--password', dest='password', help='Authentication password. Will prompt if not set.')
    parser.add_argument('--service', dest='serviceUrl', default=default_service_url,
                        help='TAP service to connect. Default: %s.' % default_service_url)
    parser.add_argument('--maxrec', dest='maxrec', default=1000000, type=int,
                        help='Maximum number of records to fecth. Default: 1000000.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='Overwrite output.')
    parser.add_argument('--force-resubmit', dest='forceResub', action='store_true',
                        help='Overwrite output.')
    parser.add_argument('--tables-dir', dest='tablesDir', default='.',
                        help='Destination path for tables. Default: current directory.')
    parser.add_argument('--job-list', dest='jobList', default='joblist.txt',
                        help='File to keel job list state. Default: ./joblist.txt')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Be verbose.')
    return parser.parse_args()

args = parse_args()

if args.debug:
    log.setLevel('DEBUG')

if args.login is None:
    login = input('Login: ')
else:
    login = args.login
if args.password is None:
    password = getpass('Password: ')
else:
    password = args.password

table_list = ['filter',
              'magabdualobj',
              'photozlephare',
              'magabsingleobj',
              'tileimage',
              'stargalclass',
              'xmatch_jplus_dr1',
              'xmatch_sdss_dr12',
              'xmatch_deep2_spec',
              #'xmatch_alhambra', # xmatch_alhambra has invalid values which mess up the VOTable parser.
              ]

log.info('Connecting to %s.' % args.serviceUrl)
tm = TAPQueueManager(args.serviceUrl, args.tablesDir)
tm.connect(login, password)

if path.exists(args.jobList):
    tm.loadJobList(args.jobList)

for tab in table_list:
    log.info('Requesting table %s.' % tab)
    tm.requestTable(tab, force=args.forceResub, maxrec=args.maxrec)
    if args.overwrite:
        tm.removeDownload(tab)

tm.saveJobList(args.jobList)

log.info('Downloading tables to %s.' % args.tablesDir)
while len(tm.listComplete()) < len(tm.jobs):
    tm.runNextJob()
    tm.downloadPending(overwrite=args.overwrite)
    log.info('Running: %s' % ','.join(tm.listRunning()))
    log.info('Still in queue: %s' % ','.join(tm.listPending()))
    sleep(20)
    
log.info('Done.')
