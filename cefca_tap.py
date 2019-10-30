'''
Created on Oct 30, 2019

@author: andre
'''

import requests
from pyvo.auth import authsession, securitymethods
import warnings
from astropy.io.votable.exceptions import VOTableSpecWarning

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


def quiet_query(service, query, maxrec=100000):
    '''
    CEFCA TAP returns VOTables that are not IVOA compliant.
    Just ignoring the silly warnings.
    Also set the maximum number of returned records to a
    more convenient value.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=VOTableSpecWarning)
        job = service.search(query, maxrec=maxrec)
    return job
