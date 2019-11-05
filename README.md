# JPAS-tables
Scripts to extract all JPAS data from VO interface to a single table.

## TAP access to CEFCA
The TAP acess provided by CEFCA requires authentication. Using `pyvo`
one needs to provide a `session` handler, which holds the cookies containing
the credentials of the connection.

The function `CEFCA_authenticate()` creates an authenticated session. It
can be used as follows:

```
from pyvo.dal import TAPService
from cefca_tap import CEFCA_authenticate
session = CEFCA_authenticate(login, password)
service = TAPService('https://archive.cefca.es/catalogues/vo/tap/minijpas-idr201910', session)
```

After that, use the `TAPService` as usual. For example, to fetch the 100 first lines
in the table `MagABDUalObj`:

```
result = service.search('select * from minijpas.MagABDualObj', maxrec=100)
tab = result.to_table()
```

The results are stored as an `astropy.table.Table` in `tab`.

Note that asynchronous jobs don't work right now, as the TAP access
does not seem to conform to the VO standards. I had to
make some changes to `pyvo` to make it work. Meanwhile,
only synchronous queries (using `TAPService.search()`as above)
are supported, which limit the maximum number of rows returned to 10000.

## Extracting tables from ADQL VO interface
When using asynchronous TAP access (by patching `pyvo`), it is
fairly simple to extract whole tables from the database. This way I download
a whole lot of tables to mess around in my computer using tools
I am already used to work. This is a poor substitution to crafting
proper ADQL scripts, but works for small amounts of data.

The asynchronous access is the same as using the VO ADQL access in
the web interface. Except that, for some reason, the list of jobs
returned by `TAPService.get_job_list()` is always empty. For one-off
jobs, one can use:

```
job = service.submit_job('select * from minijpas.MagABDualObj', maxrec=100000)
job.run()
print(job.phase)
```

After waiting the job to finish (more on this later), download the
table and save it.

```
result = job.fetch_result()
tab = result.to_table()
tab.write('my_table.fits')
```

Currently, `job.wait()`, which should block until the job is done,
is broken for this service. So, you either implement a loop, or
peek the web interface to see when it's done.

The script `download_tables.py` is a fancy version of
this procedure.

## Creating the master table
There's another script that takes the downloades tables and
produces a big table with almost all the data present.
I made this to allow an easy way to select subsamples without
having to create complicated ADQL scripts every time. The script is
`create_master_table.py`. Be aware, though, that this takes a very long time,
(a couple of hours at least) mostly to organize the single detection
magnitudes in arrays like the dual detection table.

