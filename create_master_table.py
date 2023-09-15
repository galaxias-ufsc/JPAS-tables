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
from jpdf_expand import decompress_pdf
import time

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
    parser = argparse.ArgumentParser(description='Create JPAS/JPLUS master table.')
    parser.add_argument('--output', dest='output',
                        help='Save master table to this file.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='Overwrite output.')
    parser.add_argument('--tables-dir', dest='tablesDir',
                        help='Destination path for tables.')
    parser.add_argument('--tile-id', dest='tileId',
                        help='tile_id of the tables.')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help='Be verbose.')
    return parser.parse_args()

jpas_pdf_dec_params = {'mu': [0.0, 1.5], 'sig': [0.0003333333333333333, 0.125],
                       'Nbasis': 20, 'Ncoef': 32001, 'Nmu': 751, 'Nsig': 28,
                       'Nv': 3, 'toler': 1e-10}
jplus_pdf_dec_params = {'mu': [0.0, 1.0], 'sig': [0.0008333333333333334, 0.083333333333333333],
                        'Nbasis': 20, 'Ncoef': 32001, 'Nmu': 201, 'Nsig': 103,
                        'Nv': 3, 'toler': 1e-10}

args = parse_args()
if args.debug:
    log.setLevel('DEBUG')
t0 = time.time()

# Unknown problem with data, received as 99.0
flag_99 = 1 << 31

tables_dir = args.tablesDir
tileimage_table = path.join(tables_dir, f'tileimage.fits')
mag_table = path.join(tables_dir, 'magabdualobj', f'magabdualobj_{args.tileId}.fits')
mag_single_table = path.join(tables_dir, 'magabsingleobj', f'magabsingleobj_{args.tileId}.fits')
stargalclass_table = path.join(tables_dir, 'stargalclass', f'stargalclass_{args.tileId}.fits')
photoz_table = path.join(tables_dir, 'photozlephare', f'photozlephare_{args.tileId}.fits')
mwextinction_table = path.join(tables_dir, 'mwextinction', f'mwextinction_{args.tileId}.fits')

xsdss_table = path.join(tables_dir, 'xmatch_sdss_dr12', f'xmatch_sdss_dr12_{args.tileId}.fits')

master_table = args.output
key = ['tile_id', 'number']

print('Reading dual detection table.')
master = Table.read(mag_table)
if args.debug:
    print('*** Debug: reading only 100 rows.')
    master = master[:100]
ID = ['%d-%d' % (tid, oid) for tid, oid in zip(master['tile_id'], master['number'])]
master.add_column(Column(ID), name='ID', index=0)
idx = master[key + ['ID']]


print('Reading single detection table.')
mag_single = Table.read(mag_single_table)
tileimage = Table.read(tileimage_table)
mag_single_tile = join(mag_single, tileimage, join_type='left', keys='tile_id', table_names=['ms', 'ti'])
mag_single_tile.sort(['number', 'filter_id_ms'])

# Per tile tables to speed up things.
tiles = np.unique(master['tile_id']).tolist()
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
magshape = (len(master), master['mag_auto'].shape[1])

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


for i in tqdm(range(len(master))):
    r = master[i]
    single_number = r['single_detect']
    args.tileId = r['tile_id']
    fill_single_data(data, i, single_cols, single_number, prematch[args.tileId])
            
for c, _ in single_cols:
    master[c + '_single'] = data[c]            


print('Reading xmatch SDSS table.')
xsdss = Table.read(xsdss_table)
has_zspec = xsdss['zsp'] > 0.0
xsdss_spec = xsdss[has_zspec]
xsdss_photo = xsdss[~has_zspec]
xsdss_spec_single = xmatch_get_closest(xsdss_spec, key, args.debug)
xsdss_photo_single = xmatch_get_closest(xsdss_photo, key, args.debug)

print('Reading Photo Z table.')
photoz = Table.read(photoz_table)
photoz_spec = photoz[photoz['photoz'] > 0.0]

pdf = decompress_pdf(photoz_spec['sparse_pdf'], jpas_pdf_dec_params)
cdf = np.cumsum(pdf, axis=1)
photoz_spec['z_cumulative_pdf'] = cdf

print('Reading star-gal classification table.')
stargalclass = Table.read(stargalclass_table)

print('Reading Galactic extinction table.')
mwextinction = Table.read(mwextinction_table)


print('Join SDSS.')
sdss_mag_cols = ['umag', 'gmag', 'rmag', 'imag', 'zmag']
jsdss_photo = join(idx, xsdss_photo_single, join_type='left', keys=key, table_names=['jpas', 'sdss'])
assert (jsdss_photo['tile_id'] == master['tile_id']).all()
assert (jsdss_photo['number'] == master['number']).all()
mag_sdss = np.vstack([jsdss_photo[c].filled(np.nan) for c in sdss_mag_cols])
mag_err_sdss = np.vstack([jsdss_photo['e_' + c].filled(np.nan) for c in sdss_mag_cols])
master['mag_sdss'] = mag_sdss.T
master['mag_err_sdss'] = mag_err_sdss.T
master['angdist_sdss'] = jsdss_photo['angdist']

jsdss_spec = join(idx, xsdss_spec_single, join_type='left', keys=key, table_names=['jpas', 'sdss'])
assert (jsdss_spec['tile_id'] == master['tile_id']).all()
assert (jsdss_spec['number'] == master['number']).all()
master['z_spec_sdss'] = jsdss_spec['zsp']
master['z_spec_err_sdss'] = jsdss_spec['e_zsp']
has_zsp = np.isfinite(jsdss_spec['zsp']).filled(False)
jsdss_spec = jsdss_spec[has_zsp]

#master['spobjid_sdss'] = jsdss['spobjid']
mag_sdss = np.vstack([jsdss_spec[c].filled(np.nan) for c in sdss_mag_cols])
mag_err_sdss = np.vstack([jsdss_spec['e_' + c].filled(np.nan) for c in sdss_mag_cols])
master['mag_sdss'][has_zsp] = mag_sdss.T
master['mag_err_sdss'][has_zsp] = mag_err_sdss.T
master['angdist_sdss'][has_zsp] = jsdss_spec['angdist']

print('Join Lephare photo-z.')
jphotoz = join(idx, photoz_spec, join_type='left', keys=key, table_names=['jpas', 'photoz'])
assert (jphotoz['tile_id'] == master['tile_id']).all()
assert (jphotoz['number'] == master['number']).all()
master['lephare_photoz'] = jphotoz['photoz']
master['lephare_photoz_err'] = jphotoz['photoz_err']
master['lephare_sparse_pdf'] = jphotoz['sparse_pdf']
master['lephare_z_cumulative_pdf'] = jphotoz['z_cumulative_pdf']
master['lephare_z_ml'] = jphotoz['z_ml']
master['lephare_z_best68_high'] = jphotoz['z_best68_high']
master['lephare_z_best68_low'] = jphotoz['z_best68_low']
master['lephare_chi_best'] = jphotoz['chi_best']
master['lephare_odds'] = jphotoz['odds']

print('Join star-gal classification.')
jstargalclass = join(idx, stargalclass, join_type='left', keys=key, table_names=['jpas', 'stargalclass'])
assert (jstargalclass['tile_id'] == master['tile_id']).all()
assert (jstargalclass['number'] == master['number']).all()
master['total_prob_star'] = jstargalclass['sglc_prob_star']

print('Join Galactic extinction.')
jmwextinction = join(idx, mwextinction, join_type='left', keys=key, table_names=['jpas', 'stargalclass'])
assert (jmwextinction['tile_id'] == master['tile_id']).all()
assert (jmwextinction['number'] == master['number']).all()
master['mw_ebv'] = jmwextinction['ebv']
master['mw_ebv_error'] = jmwextinction['ebv']


print(f'Writing master table: {master_table}')
master.sort(keys=['ID'])
master.write(master_table, overwrite=True)

dt = time.time() - dt
print(f'Total time: {dt:d} seconds')
