from typing import Union

import numpy as np
from astropy.io import fits
import os
from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model
from sklearn.cluster import MeanShift

dir_root = os.path.abspath(os.getcwd()+'/..')

dir_data = dir_root + '/data'
dir_code = dir_root + '/code'

def mag_to_flux(mag, magerr):
    flux = 3631e6 * 10**(-0.4*mag)
    fluxerr = flux * 0.4 * np.log(10) * magerr

    return (flux, fluxerr)

def phot_selection(phot_data):
    time, mag, magerr, flux, fluxerr = phot_data
    mask = ((magerr<0.2) & (mag<23))

    flux_data = np.array([time, flux, fluxerr])
    phot_data_good = flux_data[:, mask]
    return phot_data_good



def generate_tempfiles(band, filename, output_dir):
    data = fits.open(filename)[1].data
    data = data[data['filter']==band]

    for objid in np.unique(data['objid']):
        print objid
        data_tmp = data[data['objid'] == objid]
        mag = data_tmp['aperMag'][:,2]
        magerr = data_tmp['aperMagErr'][:,2]
        flux = data_tmp['aperFlux'][:,2]
        fluxerr = data_tmp['aperFluxErr'][:,2]
        mjd = data_tmp['mjd']

        orig_array = np.array([mjd, mag, magerr, flux, fluxerr])
        mjd, flux, fluxerr = phot_selection(orig_array)

        if len(mjd)<10:
            continue

        ### Grouping ###

        ms = MeanShift(bandwidth=0.2, bin_seeding=True)
        ms.fit(mjd.reshape(-1,1))

        labels = ms.labels_

        labels_unique = np.unique(labels)

        mask = [labels == label for label in labels_unique]
        new_mjd = [np.mean(mjd[onemask]) for onemask in mask]
        new_flux = [np.mean(flux[onemask]) for onemask in mask]
        new_fluxerr = [np.sqrt(np.sum(fluxerr[onemask]**2))/len(mjd[onemask]) for onemask in mask]

        new_array = np.array([new_mjd, new_flux, new_fluxerr]).T

        ### Output ###

        np.savetxt(output_dir + '/obj' + str(objid) + '_' +str(band), new_array, delimiter = ' ', fmt = '%.4f')




def main():

    step = 'save_lc'

    if step == 'show_light_curve':
        cfht_data = fits.open(dir_data + '/rm_phot_data/cfhtrmphot_rmqso.fits')[1].data

    elif step == 'test_javelin':
        c = get_data([dir_code + '/javelintmp/2'])
        c.plot()
        cmod = Cont_Model(c)
        cmod.do_mcmc()
        cmod.show_hist()

    elif step == 'save_lc':
        filename = dir_data + '/rm_phot_data/cfhtrmphot_rmqso.fits'
        outdir = dir_code + '/javelintmp'

        for band in 'gi':
            generate_tempfiles(band, filename, outdir)

    elif step == 'test_other':
        dummy = 1


if __name__ == '__main__':
    main()