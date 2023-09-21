'''
Created on Oct 30, 2019

@author: andre
'''

import requests
import warnings
from os import path
from astropy import log
from astropy.table import Table
import numpy as np
from time import sleep

from pyvo.auth import authsession, securitymethods
from pyvo.dal import TAPService, AsyncTAPJob

__all__ = ['CEFCA_authenticate', 'TAPQueueManager', 'download_table', 'wait_job']


def CEFCA_authenticate(login, password, login_url='https://archive.cefca.es/catalogues/login'):
    '''
    Create an authentication session to pass to ``TAPService``.
    
    Example
    -------
    
    auth = CEFCA_authenticate('my@email.com', 'p4$$w0rd')
    service = TAPService('https://archive.cefca.es/catalogues/vo/tap/minijpas-idr201910', auth)
    
    result = service.search('select * from minijpas.filter')
    print(result)
    '''
    postdata = {'login': login,
                'password': password,
                'submit': 'Sign+in',
                }

    headers = {'Content-Type': 'application/x-www-form-urlencoded',
               'Accept': 'text/plain',
               }

    session = requests.Session()
    resp = session.post(login_url, data=postdata, headers=headers)
    resp.raise_for_status()
    
    auth = authsession.AuthSession()
    auth.credentials.set(securitymethods.ANONYMOUS, session)
    return auth


def suppress_spec_warnings(func):
    from functools import wraps
    from astropy.io.votable.exceptions import VOTableSpecWarning
    from pyvo.utils.xml.exceptions import UnknownElementWarning
    _warns = [VOTableSpecWarning, UnknownElementWarning]

    @wraps(func)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            for w in _warns:
                warnings.simplefilter('ignore', category=w)
            return func(*args, **kwargs)
    return wrapper


def retry_if_error(func):
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        N = 10
        n = 1
        while n <= N:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                log.warn(f'Got exception {e}, attempt {n}/{N}.')
                n += 1
        log.error('Gave up')
        raise Exception('Too many tries.')

    return wrapper


class TAPQueueManager(object):
    '''
    Submit, run, review and fetch results of TAP jobs.
    This is a little overcomplicated because, for some reason, the job list
    returned by JPAS TAPService is always empty.
    Here we keep the job list in a local file.
    
    Parameters
    ----------
    
        service_url : str
            TAP service to connect.
            Example: https://archive.cefca.es/catalogues/vo/tap/minijpas-idr201910
        
        tables_dir : str, optional
            Directory to save tables if no path is specified when downloading.
            Default: current directory.
    '''
    
    def __init__(self, service_url, schema):
        self.serviceUrl = service_url
        self.schema = schema
        self.auth = None
        self.service = None
        self.jobs = None
        self.clearJobs()
    
    
    def connect(self, login, password):
        '''
        Authenticate the service through CEFCA portal.
        
        Parameters
        ----------
        
            login : str
                User login, usually an email.
            
            password : str
                User password.
        '''
        self.auth = CEFCA_authenticate(login, password)
        self.service = TAPService(self.serviceUrl, self.auth)
    

    def deleteAllServerJobs(self):
        log.warning('Deleting all server jobs.')
        for jd in self.service.get_job_list():
            log.warning(f'Deleting job {jd.jobid}')
            try:
                j = self._recoverJob(jd.jobid)
                j.delete()
            except Exception as e:
                log.info(f'Tried to recover and delete job {jd.jobid}. Got an exception: {e}')
    

    def _query(self, table, condition=None):
        if condition is None:
            return 'select * from %s.%s' % (self.schema, table)
        else:
            return 'select * from %s.%s where %s' % (self.schema, table, condition)
    
    
    @suppress_spec_warnings
    def syncDownloadTable(self, table, tablefile, maxrec=1000000, overwrite=False):
        '''        
        Parameters
        ----------
        
            table : str
                Name of the table, must be present in the database.
            
            tablename : str, optional
                Path to the save the table.
                
            overwrite : bool, optional
                Overwrite the file if it already exists.
                Will skip the download otherwise.
                Default: `False`

        Returns
        -------

            table_contents : Table
                Astropy Table containing the result.
        '''
        if path.exists(tablefile) and not overwrite:
            log.info(f'Table {table} already downloaded to {tablefile}, skipping download.')
            t = Table.read(tablefile)
            return t

        query = self._query(table)
        res = self.service.run_async(query, maxrec=maxrec)
        t = res_to_table(res)
        log.info('Saving table to %s.' % tablefile)
        t.write(tablefile, overwrite=overwrite)
        return t

    
    @suppress_spec_warnings
    def requestTable(self, table, tablefile, maxrec=1000000, filter=None):
        '''        
        Parameters
        ----------
        
            table : str
                Name of the table, must be present in the database.
            
            tablefile : str, optional
                Path to the save the table.
                
            maxrec : int, optional
                Maximum number of records to return.
                Default: 1000000

            filter : str, optional
                Filter to apply to query (using WHERE clause).
                Example: filter='tile_id = 123456'
                Default: None

        Returns
        -------

            job : AsyncTAPJob
                Asynchronous job queued at the server.
        '''
        query = self._query(table, filter)
        log.debug(f'Submitting query: {query}')
        j = self.service.submit_job(query, maxrec=maxrec)
        self.jobs[tablefile] = j
        log.debug(f'Starting job {j.job_id} ({table})')
        j.run()
        j.raise_if_error()
        return j
    
    
    @suppress_spec_warnings
    def download(self, tablefile, overwrite=False):
        '''
        Download results of job. The job must be already completed.
        The table will be saved as an `astropy.table.Table`, in FITS format.
        
        Parameters
        ---------
            tablefile : str
                Path to the save the table.
            
            overwrite : bool, optional
                Overwrite the file if it already exists.
                Will skip the download otherwise.
                Default: `False`
                
        '''
        if path.exists(tablefile) and not overwrite:
            log.debug('File %s already exists, skipping download.' % tablefile)
            return
        j = self.jobs[tablefile]
        log.debug('Fetching job %s results...' % j.job_id)
        res = j.fetch_result()
        
        t = res_to_table(res)
        log.info('Saving table to %s.' % tablefile)
        t.write(tablefile, overwrite=overwrite)
        log.info(f'Deleting job {j.job_id}.')
        j.delete()
        self.jobs[tablefile] = None

    
    @suppress_spec_warnings    
    def downloadPending(self, overwrite=False):
        '''
        Download all completed jobs to their default paths.
        
        Parameters
        ----------
        
            overwrite : bool, optional
                Overwrite the files if they already exist.
                Will skip the download otherwise.
                Default: `False`
            
        '''
        for tab, j in self.jobs.items():
            if j is None:
                continue
            try:
                phase = j.phase
            except Exception as e:
                log.debug(f'Tried to check job {j.job_id} phase, but got an exception.')
                log.debug(f'Exception: {e}')
                continue
            if phase == 'COMPLETED':
                log.info('Downloading table %s.' % tab)
                self.download(tab, overwrite=overwrite)
    
    
    @suppress_spec_warnings    
    def _recoverJob(self, job_id):
        url = '%s/async/%s' % (self.serviceUrl, job_id)
        job = AsyncTAPJob(url, self.auth)
        return job
        
    
    def countActiveJobs(self):
        active = [tab for tab, j in self.jobs.items() if j is not None]
        return len(active)

    
    def clearJobs(self):
        self.jobs = {}


def fix_names(t):
    '''
    Make the table columns names simpler, and transform length=1 columns to scalars.
    '''
    for name in t.colnames:
        parts = name.split('.')
        if len(parts) > 1:
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

    if isinstance(x, bytes) or isinstance(x, str):
        dtype = np.dtype(('U', max_len_x))
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
    else:
        raise Exception(f'Column type not supported: {c.name} {type(x)}.')
    return data


def res_to_table(res):
    log.debug('Converting results to table.')
    t = res.to_table()
    log.debug('Fixing table structure.')
    fix_names(t)
    convert_dtype(t)
    del t.meta['description']
    return t
