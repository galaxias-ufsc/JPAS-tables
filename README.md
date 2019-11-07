# JPAS-tables
Scripts access JPAS data from the VO TAP interface.

## TAP access to JPAS database
The TAP acess provided by CEFCA requires authentication. Using `pyvo`
one needs to provide a `session` handler, which holds the cookies containing
the credentials of the connection.

The function `CEFCA_authenticate()` creates an authenticated session. It
can be used as follows:

```
from pyvo.dal import TAPService
from jpas_tap import CEFCA_authenticate
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
This table has some issues:

* The column names all begin with `'minijpas.'`.
* All scalar columns have shape (1,).
* All array columns have dtype `Object`.

To fix these, there are some utility functions: `fix_names()` and `convert_dtype()`.
To make things easier, I created a function `download_table()` that wraps everything.
The example below downloads the table `minijpas.Filter`.

```
from jpas_tap import download_table
t = download_table(service, 'filter')
```
It may be interesting to repurpose this function for generic ADQL queries.

## Extracting tables from ADQL VO interface
When using asynchronous TAP access, it is
fairly simple to extract whole tables from the database. This way I download
a whole lot of tables to mess around in my computer using tools
I am already used to work. This is a poor substitution to crafting
proper ADQL scripts, but works fine for small amounts of data.

The asynchronous access is the same as using the VO ADQL access in
the web interface. Except that, for some reason, the list of jobs
returned by `TAPService.get_job_list()` is always empty. For one-off
jobs, one can use:

```
job = service.submit_job('select * from minijpas.MagABDualObj', maxrec=100000)
job.run()
print(job.phase)
```

To wait for the job to complete, you either do a `job.wait()`,
or implement a loop to do things while waiting (this is the whole reason
this is called asynchronous, you know.)

To download the results, do the following.

```
from jpas_tap import fix_names, convert_dtype

result = job.fetch_result()
tab = result.to_table()
fix_names(tab)
convert_dtype(tab)
tab.write('my_table.fits')
```

The script `download_tables.py` included is a fancy version of
this procedure.

## Creating the master table
There's another script that takes the downloades tables and
produces a big table with almost all the data present.
I made this to allow an easy way to select subsamples without
having to create complicated ADQL scripts every time. The script is
`create_master_table.py`. Be aware, though, that this takes a very long time,
(a couple of hours at least) mostly to organize the single detection
magnitudes in arrays like the dual detection table.

