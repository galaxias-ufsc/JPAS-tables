#!/usr/bin/env python3

'''
Created on 25 de abr de 2019

@author: andre
'''

from astropy.table import Table, Column, join
import numpy as np
from tqdm import tqdm
import argparse
from os import path
from astropy import log


def xmatch_get_closest(t, keys, debug=False):
    if debug:
        print('*** Debug: xmatch_get_closest - only 25%')
        N = len(t) // 4
        t = t[:N]
    tg = t.group_by(keys)
    rows = []
    for g in tqdm(tg.groups):
        if len(g) > 1:
            g.sort('angdist')
        rows.append(g[0].as_void())
    rows = np.array(rows, dtype=t.dtype)
    return Table(rows)

def parse_args():
    parser = argparse.ArgumentParser(description='Create JPAS master table.')
    parser.add_argument('--output', dest='output', default='master_table.fits',
                        help='Save master table to this file. Default: master-table.fits.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='Overwrite output.')
    parser.add_argument('--tables-dir', dest='tablesDir', default='.',
                        help='Destination path for tables. Default: current directory.')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Be verbose.')
    return parser.parse_args()

args = parse_args()
if args.debug:
    log.setLevel('DEBUG')

# Unknown problem with data, received as 99.0
flag_99 = 1 << 31

tables_dir = 'tables'
mag_table = path.join(tables_dir, 'magabdualobj.fits')
mag_single_table = path.join(tables_dir, 'magabsingleobj.fits')
tileimage_table = path.join(tables_dir, 'tileimage.fits')
xjplus_table = path.join(tables_dir, 'xmatch_jplus_dr1.fits')
xsdss_table = path.join(tables_dir, 'xmatch_sdss_dr12.fits')
xdeep2_table = path.join(tables_dir, 'xmatch_deep2_spec.fits')
photoz_table = path.join(tables_dir, 'photozlephare.fits')
xalhambra_table = path.join(tables_dir, 'xmatch_alhambra.fits')
muffit_table = path.join(tables_dir, 'photoz_MUFFIT_AUTO.txt')

master_table = path.join(tables_dir, 'master_v01.fits')
master_alhambra_table = path.join(tables_dir, 'master_v01_alhambra.fits')

key = ['tile_id', 'number']
jplus_fix_names= ['angdist',
                  'morph_prior_star',
                  'morph_lhood_star',
                  'morph_prob_star',
                  'gaia_prior_star',
                  'total_prob_star',
                  'lephare_photoz',
                  'lephare_photoz_err',
                  'lephare_z_cumulative_pdf',
                  'z_ml',
                  'z_best68_high',
                  'z_best68_low',
                  'chi_best',
                  'lephare_odds',
                  'tpz_photoz',
                  'tpz_photoz_err',
                  'tpz_z_cumulative_pdf',
                  'photoz_mode',
                  'photoz_mode_err',
                  'tpz_odds',
                  ]

print('Reading JPAS dual detection table.')
mag = Table.read(mag_table)
if args.debug:
    print('*** Debug: reading only 100 rows.')
    mag = mag[:100]
ID = ['%d-%d' % (tid, oid) for tid, oid in zip(mag['tile_id'], mag['number'])]
mag.add_column(Column(ID), name='ID', index=0)
idx = mag[key + ['ID']]


print('Reading JPAS single detection table.')
mag_single = Table.read(mag_single_table)
tileimage = Table.read(tileimage_table)
mag_single_tile = join(mag_single, tileimage, join_type='left', keys='tile_id', table_names=['ms', 'ti'])
mag_single_tile.sort(['number', 'filter_id_ms'])

# Per tile tables to speed up things.
tiles = np.unique(mag['tile_id']).tolist()
prematch = {}
for t in tiles:
    sel = mag_single_tile['ref_tile_id'] == t
    prematch[t] = mag_single_tile[sel]
    prematch[t].add_index('filter_id_ms')
    prematch[t].add_index('number')

single_cols = [('flags', 'int32'),
               ('mask_flags', 'int32'),
               ('fwhm_world', 'float32'),
               ('mag_auto', 'float32'),
               ('mag_err_auto', 'float32'),
               ('mag_iso', 'float32'),
               ('mag_err_iso', 'float32'),
               ('mag_petro', 'float32'),
               ('mag_err_petro', 'float32'),
               ('mu_max', 'float32'),
               ('background', 'float32'),
               ('threshold', 'float32'),
               ('mag_aper_0_8', 'float32'),
               ('mag_err_aper_0_8', 'float32'),
               ('mag_aper_1_0', 'float32'),
               ('mag_err_aper_1_0', 'float32'),
               ('mag_aper_1_2', 'float32'),
               ('mag_err_aper_1_2', 'float32'),
               ('mag_aper_1_5', 'float32'),
               ('mag_err_aper_1_5', 'float32'),
               ('mag_aper_2_0', 'float32'),
               ('mag_err_aper_2_0', 'float32'),
               ('mag_aper_3_0', 'float32'),
               ('mag_err_aper_3_0', 'float32'),
               ('mag_aper_4_0', 'float32'),
               ('mag_err_aper_4_0', 'float32'),
               ('mag_aper_6_0', 'float32'),
               ('mag_err_aper_6_0', 'float32'),
               ]
magshape = (len(mag), mag['mag_auto'].shape[1])

data = {}
for c, dt in single_cols:
    if dt == 'int32':
        fill_val = flag_99
    else:
        fill_val = np.nan
    data[c] = np.full(magshape, fill_val, dtype=dt)

def fill_single_data(data, i, single_cols, single_number, prematch):
    for j, n in enumerate(single_number):
        if n == 0: continue
        fid = j + 1
        k = np.where((prematch['filter_id_ms'] == fid) & (prematch['number'] == n))[0]
        if len(k) == 0:
            continue
        rs = prematch[k]
        if len(k) > 1:
            print('           ##### too big', k)
            rs = rs[0]
        for c, _ in single_cols:
            data[c][i, j] = rs[c]


for i in tqdm(range(len(mag))):
    r = mag[i]
    single_number = r['single_detect']
    tile_id = r['tile_id']
    fill_single_data(data, i, single_cols, single_number, prematch[tile_id])
            
for c, _ in single_cols:
    mag[c + '_jpas_single'] = data[c]            


print('Reading xmatch JPLUS table.')
xjplus = Table.read(xjplus_table)
xjplus.rename_column('jpas_tile_id', 'tile_id')
xjplus.rename_column('jpas_number', 'number')
xjplus.rename_column('jplus_tile_id', 'tile_id_jplus')
xjplus.rename_column('jplus_number', 'number_jplus')
xjplus_single = xmatch_get_closest(xjplus, key, args.debug)

print('Reading xmatch SDSS table.')
xsdss = Table.read(xsdss_table)
has_zspec = xsdss['zsp'] > 0.0
xsdss_spec = xsdss[has_zspec]
xsdss_photo = xsdss[~has_zspec]
xsdss_spec_single = xmatch_get_closest(xsdss_spec, key, args.debug)
xsdss_photo_single = xmatch_get_closest(xsdss_photo, key, args.debug)

print('Reading xmatch DEEP2 table.')
xdeep2 = Table.read(xdeep2_table)
xdeep2_spec = xdeep2[xdeep2['z'] > 0.0]
xdeep2_single = xmatch_get_closest(xdeep2_spec, key, args.debug)

print('Reading Photo Z table.')
photoz = Table.read(photoz_table)
photoz_spec = photoz[photoz['photoz'] > 0.0]

print('Join JPLUS.')
master = join(mag, xjplus_single, join_type='left', keys=key, table_names=['jpas', 'jplus'])
for c in jplus_fix_names:
    master.rename_column(c, '%s_jplus' % c)

print('Join SDSS.')
mag_cols = ['umag', 'gmag', 'rmag', 'imag', 'zmag']
jsdss_photo = join(idx, xsdss_photo_single, join_type='left', keys=key, table_names=['jpas', 'sdss'])
mag_sdss = np.vstack([jsdss_photo[c].filled(np.nan) for c in mag_cols])
mag_err_sdss = np.vstack([jsdss_photo['e_' + c].filled(np.nan) for c in mag_cols])
master['mag_sdss'] = mag_sdss.T
master['mag_err_sdss'] = mag_err_sdss.T
master['angdist_sdss'] = jsdss_photo['angdist']

jsdss_spec = join(idx, xsdss_spec_single, join_type='left', keys=key, table_names=['jpas', 'sdss'])
master['z_spec_sdss'] = jsdss_spec['zsp']
master['z_spec_err_sdss'] = jsdss_spec['e_zsp']
has_zsp = np.isfinite(jsdss_spec['zsp']).filled(False)
jsdss_spec = jsdss_spec[has_zsp]

#master['spobjid_sdss'] = jsdss['spobjid']
mag_sdss = np.vstack([jsdss_spec[c].filled(np.nan) for c in mag_cols])
mag_err_sdss = np.vstack([jsdss_spec['e_' + c].filled(np.nan) for c in mag_cols])
master['mag_sdss'][has_zsp] = mag_sdss.T
master['mag_err_sdss'][has_zsp] = mag_err_sdss.T
master['angdist_sdss'][has_zsp] = jsdss_spec['angdist']


print('Join DEEP2 spectra.')
jdeep2 = join(idx, xdeep2_single, join_type='left', keys=key, table_names=['jpas', 'deep2'])
master['angdist_deep2'] = jdeep2['angdist']
master['z_spec_deep2'] = jdeep2['z']
master['z_spec_err_deep2'] = jdeep2['e_z']

print('Join Lephare photo-z.')
jphotoz = join(idx, photoz_spec, join_type='left', keys=key, table_names=['jpas', 'photoz'])
master['z_photo_lephare'] = jphotoz['photoz']
master['z_photo_err_lephare'] = jphotoz['photoz_err']

if path.exists(muffit_table):
    print('Reading MUFFIT Photo Z table.')
    muffit = Table.read(muffit_table, format='ascii', names=['ID', 'ra', 'dec', 'z_photo', 'z_photo_err'])
    
    print('Join MUFFIT.')
    jmuffit = join(idx, muffit, join_type='left', keys='ID', table_names=['jpas', 'muffit'])
    master['z_photo_muffit'] = jmuffit['z_photo']
    master['z_photo_err_muffit'] = jmuffit['z_photo_err']

print('Writing master table.')
master.sort(keys=['ID'])
master.write(master_table, overwrite=True)

if path.exists(xalhambra_table):
    print('Reading Alhambra.')
    xalhambra = Table.read(xalhambra_table)
    
    alhambra_filters = ['f396w',
                        'f427w',
                        'f458w',
                        'f489w',
                        'f520w',
                        'f551w',
                        'f582w',
                        'f613w',
                        'f644w',
                        'f675w',
                        'f706w',
                        'f737w',
                        'f768w',
                        'f799w',
                        'f830w',
                        'f861w',
                        'f892w',
                        'f923w',
                        'f954w',
                        ]
    
    mag_shape = (len(xalhambra), len(alhambra_filters))
    alhambra_mag = np.zeros(mag_shape)
    alhambra_mag_err = np.zeros(mag_shape)
    alhambra_irms = np.zeros(mag_shape)
    
    for i, f in enumerate(alhambra_filters):
        alhambra_mag[:, i] = xalhambra[f]
        alhambra_mag_err[:, i] = xalhambra['e_' + f]
        alhambra_irms[:, i] = xalhambra['irms' + f]
    
        xalhambra.remove_column(f)
        xalhambra.remove_column('e_' + f)
        xalhambra.remove_column('irms' + f)
        xalhambra.remove_column('l_' + f)
    
    xalhambra['mag'] = alhambra_mag
    xalhambra['mag_err'] = alhambra_mag
    xalhambra['irms'] = alhambra_irms
    
    for c in xalhambra.colnames:
        if c in ['tile_id', 'number', 'angdist']: continue
        xalhambra.rename_column(c, '%s_alhambra' % c)
    
    xalhambra_single = xmatch_get_closest(xalhambra, key, args.debug)
    
    print('Join Alhambra.')
    master_alhambra = join(master, xalhambra_single, join_type='left', keys=key, table_names=['jpas', 'alhambra'])
    
    
    print('Writing master Alhambra table.')
    master_alhambra.write(master_alhambra_table, overwrite=True)
