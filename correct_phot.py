import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os, sys

# import astroML as ml

dir_root = os.getcwd() + '/../..'
dir_code = dir_root + '/code'
dir_data = dir_root + '/data'


def check_cadence():
    filename = dir_data + '/SDSS-RM/rm_phot_data/cfhtrmphot_rmqso.fits'
    data = fits.open(filename)[1].data
    print(len(data))
    data = data[data['aperMagErr'][:, 1] < 0.05]

    # print(data.columns)

    for rmid in range(10, 20):
        data_thisq = data[data['objId'] == rmid]
        print(len(data_thisq))
        plt.errorbar(data_thisq['mjd'], data_thisq['aperMag'][:, 1], yerr=data_thisq['aperMagErr'][:, 1], fmt='.')
        plt.show()


check_cadence()
