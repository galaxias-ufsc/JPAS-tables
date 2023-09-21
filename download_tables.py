'''
Created on Oct 30, 2019

@author: andre
'''

import argparse
from getpass import getpass
from jpas_tap import TAPQueueManager
from astropy import log
from astropy.table import Table
from os import path, makedirs
import time

default_service_url = 'https://archive.cefca.es/catalogues/vo/tap/minijpas-pdr201912'
default_schema = 'minijpas'

def parse_args():
    parser = argparse.ArgumentParser(description='Download JPAS tables by tile.')
    parser.add_argument('tileList', type=str, nargs='*',
                        help='Download tables for these tiles.')
    parser.add_argument('--login', dest='login', help='Authentication username. Will prompt if not set.')
    parser.add_argument('--password', dest='password', help='Authentication password. Will prompt if not set.')
    parser.add_argument('--service', dest='serviceUrl', default=default_service_url,
                        help='TAP service to connect. Default: %s.' % default_service_url)
    parser.add_argument('--schema', dest='schema', default=default_schema,
                        help='Database schema. Default: %s.' % default_schema)
    parser.add_argument('--maxrec', dest='maxrec', default=1000000, type=int,
                        help='Maximum number of records to fecth. Default: 1000000.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='Overwrite output.')
    parser.add_argument('--tables-dir', dest='tablesDir', default='.',
                        help='Destination path for tables. Default: current directory.')
    parser.add_argument('--job-list', dest='jobList', default='joblist.txt',
                        help='File to keel job list state. Default: ./joblist.txt')
    parser.add_argument('--delete-all-server-jobs', dest='deleteAllServerJobs', action='store_true',
                        help='Delete all server jobs on start.')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Be verbose.')
    return parser.parse_args()

def get_ref_tiles(tileimage_file):
    t = Table.read(tileimage_file)
    return list(set(t['ref_tile_id']))


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


meta_table_list = ['tileimage',
                   'filter',
                   ]

tileimage = 'tileimage'


big_table_list = ['magabdualobj',
                 'photozlephare',
                 'stargalclass',
                 'mwextinction',
                 'xmatch_sdss_dr12',
                 ]

log.info('Connecting to %s.' % args.serviceUrl)
tm = TAPQueueManager(args.serviceUrl, args.schema)
tm.connect(login, password)

if args.deleteAllServerJobs:
    tm.deleteAllServerJobs()

log.debug('Creating destination directories.')
if not path.exists(args.tablesDir):
    makedirs(args.tablesDir)
for table in big_table_list:
    td = path.join(args.tablesDir, table)
    if not path.exists(td):
        makedirs(td)

log.info('Downloading metadata tables.')
for table in meta_table_list:
    tablefile = path.join(args.tablesDir, f'{table}.fits')
    log.info(f'Downloading table {table} to {tablefile}.')
    tm.syncDownloadTable(table, tablefile, overwrite=args.overwrite)

if len(args.tileList) > 0:
    log.info('Using tiles listed in command line.')
    ref_tiles = args.tileList
else:
    log.info('No tiles listed in command line.')
    tileimage_file = path.join(args.tablesDir, f'{tileimage}.fits')
    log.info(f'Downloading all tiles from {tileimage_file}.')
    ref_tiles = get_ref_tiles(tileimage_file)
tile_count = 1
for tid in ref_tiles:
    t0 = time.time()
    log.info(f'Current tile_id: {tid} ({tile_count}/{len(ref_tiles)})')
    for table in big_table_list:
        tablefile = path.join(args.tablesDir, table, f'{table}_{tid}.fits')
        if path.exists(tablefile) and not args.overwrite:
            log.info(f'Table {tablefile} already downloaded.')
            continue
        log.info(f'Requesting table {table}.')
        j = tm.requestTable(table, tablefile, filter=f'tile_id={tid}', maxrec=args.maxrec)
        log.debug(f'Created job {j.job_id} for table {tablefile}.')
    while tm.countActiveJobs() > 0:
        time.sleep(10)
        tm.downloadPending(overwrite=args.overwrite)
        log.info(f'{tm.countActiveJobs()} jobs still in queue')
        log.debug('Jobs in queue: %s' % ','.join(tm.jobs.keys()))
    tile_count += 1
    dt = time.time() - t0
    log.info(f'Took {dt:.0f} seconds to download data for tile {tid}.')

log.info('Done.')
