### Correct all filters to CFHT g or i ####

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

from math import *
from astropy.io import fits
from scipy import integrate
from astropy import units as u
from astropy import constants as const
from sklearn.cluster import MeanShift


dir_root = os.path.abspath(os.getcwd()+'/..')

dir_data = dir_root + '/data'
dir_code = dir_root + '/code'

data_file_dir = {'cfht': dir_data + '/rm_phot_data/cfhtrmphot_allqso.fits',\
                 'bok': dir_data + '/rm_phot_data/bokrmphot_allqso.fits',\
                 'sdss': dir_data + '/rm_phot_data/allqso_rmfield_xmatched.fits',\
                 'pstar': dir_data + '/rm_phot_data/panstarrs_stackobject_search.txt'}


class lambdafunc(object):
    def update(self, wavelength, value, units):
        if not len(wavelength) == len(value):
            raise ValueError('The length of wavelength list does not equal to the length of value list')
        else:
            self.wavelength = np.array(wavelength) * float(1 * units[0].cgs / u.Angstrom)
            self.value = np.array(value)
            self.units = [u.Angstrom.cgs, units[1].cgs]
            self.interval = (np.append(self.wavelength, 0) - np.append(0, self.wavelength))[1:-1]
            self.midvalue = (np.append(self.value, 0) + np.append(0, self.value))[1:-1] / 2.0

    def getvalue(self, x):
        return np.interp(x, self.wavelength, self.value)

    def plot(self):
        plt.plot(self.wavelength, self.value)
        plt.xlabel('Wavelength / $\AA$', fontsize=16)
        plt.ylabel('Arbitrary Unit', fontsize=16)
        plt.show()


class FiltCurve(lambdafunc):
    def __init__(self, wavelength, value, units=[u.Angstrom, u.Quantity(1)]):
        self.update(wavelength, value, units)


class Spectrum(lambdafunc):
    def __init__(self, wavelength, value, units=[u.Angstrom, u.erg / u.s / u.cm / u.cm / u.Angstrom], mode='OBS'):
        self.update(wavelength, value, units)
        self.mode = mode

    def flux(self, filt):
        newflux = self.getvalue(np.array(filt.wavelength * filt.units[0] / self.units[0]))
        filtered_flux = newflux * filt.value
        objflux = integrate.simps(filtered_flux, x=filt.wavelength)

        return objflux


    def magnitude(self, filt, style='AB'):
        '''
        Calculate the magnitude in the input filter.
        '''
        ### Now, AB mag only ###

        objflux = self.flux(filt)

        stdspec = ((3631 * u.Jy * const.c).cgs / ((filt.wavelength * filt.units[0]).cgs) ** 2).to(self.units[1])
        stdflux = integrate.simps(stdspec, x=filt.wavelength)

        mag = -2.5 * np.log10(float(objflux / stdflux))
        return mag


    def normalize(self, filt, mag, style='AB'):
        mag0 = self.magnitude(filt, style)
        scale = 10 ** (0.4 * (mag0 - mag))

        self.value = self.value * scale



class lightcurve(object):
    def __init__(self, sdss_file, ps_file, bok_file, cfht_file, spec_file, info_file):
        dummy = 1


def grouping_mjd(mjd):
    ms = MeanShift(bandwidth=0.2, bin_seeding=True)
    ms.fit(mjd.reshape(-1, 1))

    labels = ms.labels_

    labels_unique = np.unique(labels)

    mask = [labels == label for label in labels_unique]
    return mask


def mag_to_flux(mag, magerr):
    flux = 3631e6 * 10**(-0.4*mag)
    fluxerr = flux * 0.4 * np.log(10) * magerr

    ### meanshift ###

    return (flux, fluxerr)


def read_psdata(filename):
    f = open(filename)
    origID = 0
    # for idx in range(3):
    line0 = f.readline()
    line1 = f.readline()

    line = f.readline()

    IDlist = []
    contentlist = []

    while len(line) > 0:
        line = f.readline()
        if len(line) == 0:
            break
        if 'Input line' in line:
            origID += 1

        else:
            IDlist.append(str(origID))
            contentlist.append(line)
        # print(line)

    np.savetxt('./pstmp', np.transpose([IDlist, contentlist]), delimiter=',', fmt='%500s')

    data = pd.read_csv('./pstmp', names=['origID'] + line1.split(','))
    os.system('rm ./pstmp')
    return data



class QSO(object):

    def __init__(self, objid, info_file = dir_data + '/rm_data_files/allqsos_rmfield.fits'):
        allinfo = fits.open(info_file)[1].data
        this_quasar_info = allinfo[allinfo['objid']==objid]

        # basic information #

        self.objid = objid
        self.redshift = this_quasar_info['Z']
        self.coord = (this_quasar_info['ra'], this_quasar_info['dec'])

        # load spectrum #

        #self.spectrum = Spectrum(wavelength, flux)


    def read_phot(self, band, key):

        #sdss_file, ps_file, bok_file, cfht_file
        #mag_key = band + 'mag'

        data_file = data_file_dir[key]

        objid = self.objid

        programs = ('cfht', 'bok', 'sdss', 'pstar')
        if not key in programs:
            raise ValueError('The program name is not recognized. Only %s, %s, %s and %s\
                             data are available' %programs)

        if key == 'cfht':

            # CHFT data #
            cfht = fits.open(data_file)[1].data
            cfht_obj = cfht[cfht['objid']==objid]
            cfht_band = cfht_obj[cfht_obj['filter']==band]

            # Use the nth aperture #
            n_aperture_cfht = 3
            cfht_mjd = cfht_band['mjd']
            cfht_mag = cfht_band['aperMag'][:, n_aperture_cfht]
            cfht_mag_error = cfht_band['aperMagErr'][:, n_aperture_cfht]

            # Meanshift #
            cfht_mask = grouping_mjd(cfht_mjd)

            cfht_new_mjd = [np.mean(cfht_mjd[onemask]) for onemask in cfht_mask]
            cfht_new_mag = [np.mean(cfht_mag[onemask]) for onemask in cfht_mask]
            cfht_new_magerr = [np.sqrt(np.sum(cfht_mag_error[onemask] ** 2)) /\
                          len(cfht_mjd[onemask]) for onemask in cfht_mask]

            self.cfht_phot = {band: {'mjd': cfht_new_mjd,\
                              'mag': cfht_new_mag,\
                              'err': cfht_new_magerr}}

        elif key == 'bok':

            # Bok Data #
            bok = fits.open(data_file)[1].data
            bok_obj = bok[bok['objid'] == objid]
            bok_band = bok_obj[bok_obj['filter']==band]

            # Use the nth aperture #
            n_aperture_bok = 3
            bok_mjd = bok_band['mjd']
            bok_mag = bok_band['aperMag'][:, n_aperture_bok]
            bok_mag_error = bok_band['aperMagErr'][:, n_aperture_bok]

            # Meanshift #
            bok_mask = grouping_mjd(bok_mjd)
            bok_new_mjd = [np.mean(bok_mjd[onemask]) for onemask in bok_mask]
            bok_new_mag = [np.mean(bok_mag[onemask]) for onemask in bok_mask]
            bok_new_magerr = [np.sqrt(np.sum(bok_mag_error[onemask] ** 2)) /\
                          len(bok_mjd[onemask]) for onemask in bok_mask]

            self.bok_phot = {band: {'mjd': bok_new_mjd, \
                              'mag': bok_new_mag, \
                              'err': bok_new_magerr}}


        elif key == 'pstar':

            # Pan-Starrs Data #
            pstar = read_psdata(data_file)
            pstar_obj = pstar[pstar['origID']==objid]
            pstar_mjd = pstar_obj['epochMean']
            pstar_mag = pstar_obj[band + 'PSFMag']
            pstar_magerr = pstar_obj[band + 'PSFMagErr']

            self.pstar_phot = {band: {'mjd': pstar_mjd, \
                              'mag': pstar_mag, \
                              'err': pstar_magerr}}

        elif key == 'sdss':

            # SDSS Data #
            sdss = fits.open(data_file)[1].data
            sdss_obj = sdss[sdss['m_OBJID']==objid]
            sdss_mjd = sdss_obj['s_mjd']
            sdss_mag = sdss_obj['s_psfMag_'+band]
            sdss_magerr = sdss_obj['s_psfMagErr_'+band]

            self.sdss_phot = {band: {'mjd': sdss_mjd, \
                              'mag': sdss_mag, \
                              'err': sdss_magerr}}


    def merge_phot(self, band):

        programs = ('cfht', 'bok', 'sdss', 'pstar')

        # Read raw photometry data #
        for key in programs:
            self.read_phot(band, key)

        # Filter correction to CFHT MegaCam g #
        # Bok #


        """
        # Filter Correction #
        specs = fits.open(spec_file)[1].data
        infos = fits.open(info_file)[1].data
    
        for objid in infos['objID']:
            spec_thisq = specs[objid]
            wavelength, flux, fluxerr = spec_thisq
            specobj_thisq = QSO_spectrum(wavelength, flux)
    
        """
def main():
    qso0 = QSO(0)
    qso0.read_phot('g', 'sdss')
    print(qso0.sdss_phot)


if __name__ == '__main__':
    main()