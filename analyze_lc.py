from typing import Union

import numpy as np
from astropy.io import fits
import os
from javelin.zylc import get_data

dir_root = os.path.abspath(os.getcwd()+'/..')

dir_data = dir_root + '/data'

def main():

    step = 'show_light_curve'

    if step == 'show_light_curve':
        cfht_data = fits.open(dir_data + '/rm_phot_data/cfhtrmphot_rmqso.fits')[1].data


if __name__ == '__main__':
    main()