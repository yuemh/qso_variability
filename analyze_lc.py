from typing import Union

import numpy as np
from astropy.io import fits
import os
from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model

dir_root = os.path.abspath(os.getcwd()+'/..')

dir_data = dir_root + '/data'

def main():

    step = 'test_javelin'

    if step == 'show_light_curve':
        cfht_data = fits.open(dir_data + '/rm_phot_data/cfhtrmphot_rmqso.fits')[1].data

    elif step == 'test_javelin':
        c = get_data(['/Users/minghao/Software/javelin-0.33/examples/dat/continuum.dat'])
        print(c)


if __name__ == '__main__':
    main()