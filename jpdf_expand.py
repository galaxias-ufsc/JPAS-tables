#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) (2020) CEFCA
#
# jpdf-expand is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# jpdf-expand is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with jpdf-expand. If not, see <http://www.gnu.org/licenses/>.
#
# This script is based on the SparsePz module created by
# Matias Carrasco Kind https://github.com/mgckind/SparsePz
#

from __future__ import print_function
import sys
from numpy import *
from scipy import linalg as sla
from scipy.optimize import leastsq
from scipy import special
from astropy.io import fits
from astropy.table import Table, Column
import copy
import argparse

def decompress_pdf(sparsePDF, params):
    """
    Extract the full PDFs from the indices of the sparse representation

    parameters:
    sparsePDF: N x M array containing the M integer indices for each of the N galaxies
    params: dictionary containing values for configuration parameters used in compression.

    returns:
    N x Nz array of floats containing the decompressed PDFs, where N is the number of
    galaxies and Nz is the number of redshift values in which each PDF is sampled.
    """
    zstep = (params['mu'][1]-params['mu'][0])/(params['Nmu']-1)
    zv = arange(params['mu'][0],params['mu'][1]+0.00001,zstep,dtype=float)
    pars = copy.copy(params)
    pars['z'] = zv
    Nsources, Nbasis = sparsePDF.shape
    outPDF = zeros((Nsources,len(zv)),dtype=float)

    for k in range(Nsources):
        sp = sparsePDF[k][:] # long indices
        while (len(sp) > 0) and (sp[-1] == 0):
            sp = sp[:-1]  # remove trailing zeros

        if len(sp) > 0:
            outPDF[k] = reconstruct_pdf_int(sp, pars, cut=1.e-5)
        else:
            pass # no data for this source

    return outPDF

def reconstruct_pdf_int(long_index, header, cut=1.e-5):
    """
    This function reconstruct the pdf from the integer indices only and the parameters used to create the dictionary
    with Gaussians and Voigt profiles

    :param int long_index: List of indices including coefficients (32bits integer array)
    :param dict header: Dictionary of the fits file header with information used to create dictionary and sparse indices
    :param float cut: cut threshold when creating the dictionary

    :return: the pdf normalized so it sums to one
    """

    Ncoef = header['Ncoef']
    zfine = header['z']
    mu = header['mu']
    Nmu = header['Nmu']
    sigma = header['sig']
    Nsigma = header['Nsig']
    Nv = header['Nv']

    VALS = linspace(0, 1, Ncoef)
    dVals = VALS[1] - VALS[0]
    sp_ind = array(list(map(get_N, long_index)))
    spi = sp_ind[:, 0]
    Dind2 = sp_ind[:, 1]
    vals = spi * dVals
    ####
    vals[0]=1.
    ####

    rep_pdf = reconstruct_pdf_v(Dind2, vals, zfine, mu, Nmu, sigma, Nsigma, Nv)

    return rep_pdf


def reconstruct_pdf_v(index, vals, zfine, mu, Nmu, sigma, Nsigma, Nv, cut=1.e-5):
    """
    This function reconstruct the pdf from the indices and values and parameters used to create the dictionary with
    Gaussians and Voigt profiles

    :param int index: List of indices in the dictionary for the selected bases
    :param float vals: values or coefficients corresponding to the listed indices
    :param float zfine: redshift values from the original pdf or used during the sparse representation
    :param float mu: [min_mu, max_mu] values used to create the dictionary
    :param int Nmu: Number of mu values used to create the dictionary
    :param float sigma: [min_sigma, mas_sigma] sigma values used to create the dictionary
    :param int Nsigma: Number of sigma values
    :param int Nv: Number of Voigt profiles used to create dictionary
    :param float cut: cut threshold when creating the dictionary

    :return: the pdf normalized so it sums to one
    """

    zmid = linspace(mu[0], mu[1], Nmu)
    sig = linspace(sigma[0], sigma[1], Nsigma)
    gamma = linspace(0, 0.5, Nv)
    pdf = zeros(len(zfine))
    for kk in range(len(index)):
        i = int(floor(index[kk] / (Nsigma * Nv)))
        j = int(floor((index[kk] % (Nsigma * Nv)) / Nv))
        k = (index[kk] % (Nsigma * Nv)) % Nv

        pdft = voigt(zfine, zmid[i], sig[j], sig[j] * gamma[k])

        pdft = where(pdft >= cut, pdft, 0.)
        pdft = pdft / linalg.norm(pdft)
        pdf += pdft * vals[kk]

    pdf = where(greater(pdf, max(pdf) * 0.005), pdf, 0.)
    if sum(pdf) > 0: pdf = pdf / sum(pdf)
    return pdf

def voigt(x, x_mean, sigma, gamma):
    """
     Voigt profile
     V(x,sig,gam) = Re(w(z)), w(z) Faddeeva function
     z = (x+j*gamma)/(sigma*sqrt(2))

     :param float x: the x-axis values (redshift)
     :param float x_mean: Mean of the gaussian or Voigt
     :param float sigma: Sigma of the original Gaussian when gamma=0
     :param float gamma: Gamma parameter for the Lorentzian profile (Voigt)

     :return: The real values of the Voigt profile at points x
     """

    x = x - x_mean
    z = (x + 1j * gamma) / (sqrt(2.) * sigma)
    It = special.wofz(z).real
    return It

def get_N(longN):
    """
    Extract coefficients fro the 32bits integer,
    Extract Ncoef and Nbase from 32 bit integer
    return (longN >> 16), longN & 0xffff

    :param int longN: input 32 bits integer

    :return: Ncoef, Nbase both 16 bits integer
    """
    return (longN >> 16), (longN & (2 ** 16 - 1))


if __name__ == '__main__':
    # configure argument parser
    parser = argparse.ArgumentParser(description='''Reconstruct the photoz probability
        distribution function (PDF) of miniJPAS sources from compressed representations.''',
        epilog='''The input table must include
        a column named "SPARSE_PDF", containing the sparse indices that enconde each
        source's full photoz PDF in the miniJPAS database.
        The output is a modified version of the input table where the column "SPARSE_PDF"
        is replace by "PDF", which contains the full (decompressed) PDF of the sources.''')

    parser.add_argument('input', metavar='INPUT.fits', type=str,
                        help='name of input FITS table containing compressed PDFs')
    parser.add_argument('output', metavar='OUT.fits', type=str, nargs='?',
    		    help='(optional) name for output FITS table. If ommited, the input table is overwritten.')

    args = parser.parse_args()

    if args.output == None:
        args.output = args.input # overwrite input table

    params = {'mu': [0.0, 1.5], 'sig': [0.0003333333333333333, 0.125],
              'Nbasis': 20, 'Ncoef': 32001, 'Nmu': 751, 'Nsig': 28,
              'Nv': 3, 'toler': 1e-10}

    hdu_list = fits.open(args.input)
    data = hdu_list[1].data
    cols_orig = hdu_list[1].columns.names
    cols = array([str(i).upper() for i in cols_orig]) # uppercase names
    table = Table(data)

    if 'SPARSE_PDF' in cols:
        sparse = asarray(data['SPARSE_PDF'])
        print("decompressing PDFs...")
        exp_PDF = decompress_pdf(sparse, params)
        print("done.")
        PDF = Column(exp_PDF, name='PDF')

        index = where(cols == 'SPARSE_PDF')[0]
        table.remove_column(cols_orig[index[0]])
        table.add_column(PDF,index=index[0]) # index indicates position where column is inserted
    else:
        print('ERROR: column "SPARSE_PDF" not found in table')
        sys.exit(1)

    print("writting to table {0}...".format(args.output))
    table.write(args.output, format='fits', overwrite='True')
    print("done.")
