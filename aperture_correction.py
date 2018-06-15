### Aperture Correction ###

import numpy as np
import os
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.integrate import quad
from math import *

def moffat_flux(r, para):
    I, alpha, beta = para
    result = quad(lambda x: I * (beta - 1)/pi/alpha**2 * (1+(x/alpha)**2)**(-beta), 0, r)

    return result.y






def fit_profile(data, model):
    radius, flux = data
    result = curve_fit(moffat_flux, radius, flux, [flux[-1], 1, 4])

    return result.x


def main():
    dummy = 1

if __name__ == '__main__':
    main()
