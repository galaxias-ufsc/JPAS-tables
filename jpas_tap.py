'''
Created on Oct 30, 2019

@author: andre
'''

import requests
import warnings
from os import path
from astropy import log
import numpy as np


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
    
    def __init__(self, service_url, tables_dir='.'):
        self.serviceUrl = service_url
        self.tablesDir = tables_dir
        self.auth = None
        self.service = None
        self.jobs = {}
    
    
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
    
    
    def loadJobList(self, filename):
        '''
        Recover job information from a list kept in an ascii table in `filename`.
        The lines must contain:

            table_name job_id PHASE
            
        where job_id is an integer, and PHASE is a (informative only) string.
        
        Parameters
        ----------
        
            filename : str
                Ascii file containing the submitted jobs.
        '''
        with open(filename, 'r') as f:
            for l in f.readlines():
                w = l.split()
                if len(w) != 3:
                    continue
                tablename, job_id, phase = w
                log.debug('Recovering table %s (job id %s, phase %s)' % (tablename, job_id, phase))
                self._recoverJob(tablename, job_id)
    
    
    @suppress_spec_warnings
    def saveJobList(self, filename):
        '''
        Save current submitted jobs information.
        NOTE: this will overwrite the file.
        
        Paramters
        ---------
            filename : str
                Ascii file containing the submitted jobs.
        
        '''
        with open(filename, 'w') as f:
            f.writelines('%s %s %s\n' % (t, j.job_id, j.phase) for t, j in self.jobs.items())
    
    
    def _query(self, tablename):
        return 'select * from minijpas.%s' % tablename
    
    
    def _filePath(self, tablename):
        return path.join(self.tablesDir, '%s.fits' % tablename)
    
    
    @suppress_spec_warnings
    def requestTable(self, tablename, force=False, maxrec=1000000):
        '''
        
        Parameters
        ----------
        
            tablename : str
                Name of the table, must be present in the database.
            
            force : bool, optional
                If the table is already in the queue, the default
                behaviour is to use that job. if `force=True`, submit another job.
            
            maxrec : int, optional
                Maximum number of records to return.
                Default: 1000000
        '''
        if tablename in self.jobs and not force:
            log.debug('Table %s already submitted. Use force=True to resubmit.' % tablename)
            return
        query = self._query(tablename)
        j = self.service.submit_job(query, maxrec=maxrec)
        log.debug('Created job %s for table %s.' % (j.job_id, tablename))
        self.jobs[tablename] = j
        return j
    
    
    @suppress_spec_warnings
    def runNextJob(self):
        '''
        Run next pending job.
        '''
        for tab, j in self.jobs.items():
            if j.phase == 'PENDING':
                j.run()
                log.info('Started job %s (%s).' % (j.job_id, tab))
                return j
        return None
    
    
    @suppress_spec_warnings
    def _listFiltered(self, phases):
        return [t for t, j in self.jobs.items() if j.phase in phases]
    
    
    def listComplete(self):
        '''
        List complete jobs, which have results ready for download.
        
        Returns
        -------

            complete_jobs : list
                List containing the table names of the complete jobs.
        '''
        return self._listFiltered(['COMPLETED'])
    
    
    def listRunning(self):
        '''
        List jobs currently running.

        Returns
        -------

            complete_jobs : list
                List containing the table names of the running jobs.
        '''
        return self._listFiltered(['RUN', 'EXECUTING'])
    
    
    def listPending(self):
        '''
        List jobs currently in queue.

        Returns
        -------

            complete_jobs : list
                List containing the table names of the queued jobs.
        '''
        return self._listFiltered(['PENDING'])
    
    
    def _downloaded(self, tablename):
        filename = self._filePath(tablename)
        return path.exists(filename)
        

    def removeDownload(self, tablename):
        filename = self._filePath(tablename)
        if path.exists(filename):
            log.debug('Deleting file %s.' % tablename)
        
        
    def listDownloadPending(self):
        '''
        List complete jobs, not yet downloaded to the default path.
        

        Returns
        -------

            available_jobs : list
                List containing the table names of the downloadable jobs.
        '''
        return [tab for tab in self.listComplete() if not self._downloaded(tab)]
    
    
    @suppress_spec_warnings
    def download(self, tablename, filename=None, overwrite=False):
        '''
        Download results of job. The job must be already completed.
        The table will be saved as an `astropy.table.Table`, in FITS format.
        
        Parameters
        ---------
            tablename : str
                Name of the table to download.
            
            filename : str, optional
                Path to the save the table.
                Default: `'table_dir/tablename.fits'`, where
                `table_dir` may be set in the constructor(falls back to `'.'`).
                
            overwrite : bool, optional
                Overwrite the file if it already exists.
                Will skip the download otherwise.
                Default: `False`
                
        '''
        if filename is None:
            filename = self._filePath(tablename)
        if path.exists(filename) and not overwrite:
            log.debug('File %s already exists, skipping download.' % filename)
            return
            
        j = self.jobs[tablename]
        if j.phase != 'COMPLETED':
            raise Exception('Job %s (table %s) is in phase %s. Download unavailable.' % (j.job_id, tablename, j.phase))
        
        log.debug('Fetching job %s results...' % j.job_id)
        res = j.fetch_result()
        
        j.raise_if_error()
        t = res.to_table()
        log.debug('Fixing table structure.')
        fix_names(t)
        convert_dtype(t)
        del t.meta['description']
        log.info('Saving table to %s.' % filename)
        t.write(filename, overwrite=overwrite)
    
    
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
        for tab in self.listDownloadPending():
            log.info('Downloading table %s.' % tab)
            self.download(tab, overwrite=overwrite)
    
    
    @suppress_spec_warnings    
    def _recoverJob(self, tablename, job_id):
        url = '%s/async/%s' % (self.serviceUrl, job_id)
        job = AsyncTAPJob(url, self.auth)
        self.jobs[tablename] = job
        

@suppress_spec_warnings    
def download_table(service, table_name, filename=None, overwrite=False, maxrec=1000000, timeout=600.0):
    '''
    Download a table from a TAP service, fix the columns and
    save to disk. Uses async mode to allow more records to be retrieved.
    
    Parameters
    ----------
        service : TAPService
        
        table_name : string
            Table name without the schema prefix. (Ex. no 'minijpas.')
        
        filename : string, optional
            If set, write the table to this path to the file to be saved.
            Must have a proper extension (.fits, etc.) to allow astropy to guess the format.
            
        overwrite : bool, optional
            Overwrite the file if it already exists.
            Default: `False`

        maxrec : int, optional
            Maximum number of records to return.
            Default: 1000000

        timeout : float, optional
            Timeout in seconds to the wait loop while waiting for the job to complete.
            Default: 600.0

    Returns
    -------
        table : astropy.table.Table
            A table containing the results.
    '''
    if path.exists(filename) and not overwrite:
        log.debug('File %s already exists, skipping download.')
        return
    query = 'select * from minijpas.%s' % table_name

    log.info('Submitting query.')
    log.debug('Query: %s' % query)
    job = service.submit_job(query, maxrec=maxrec)
    log.info('Starting job.')
    job.run()
    log.info('Waiting the job completion (current phase: %s)...' % job.phase)
    job.wait(timeout=timeout)
    log.info('Fetching job results...')
    result = job.fetch_result()
        
    job.raise_if_error()
    t = result.to_table()
    log.debug('Fixing table structure.')
    fix_names(t)
    convert_dtype(t)
    del t.meta['description']
    if filename is not None:
        log.info('Saving table to %s.' % filename)
        t.write(filename, overwrite=overwrite)
    return t


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