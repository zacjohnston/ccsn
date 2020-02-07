import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# adapted from https://github.com/snaphu-msu/ecRateStudy


def get_bounce_time(masses, filenames):
    bouncetimes = {}

    for mass in masses:
        bouncetimes[mass] = {}
        logname = filenames[mass][:-4] + ".log"

        for line in open(logname, 'r'):
            if "Bounce!" in line:
                bouncetimes[mass]["tbounce"] = float(line.split()[1])

        print(mass, bouncetimes[mass]["tbounce"])

    return bouncetimes


def extract_last_dats(masses, filenames, var_list=None):
    """Extract last line of .dat files
    """
    print('Extracting last lines of /dat files')

    if var_list is None:
        var_list = {'time': 0,
                    'exp_en': 9,
                    'rsh_avg': 11,
                    'dens_c': 16,
                    'pns_mass': 20}

    last_dats = pd.DataFrame()
    last_dats['mass'] = masses

    # create temporary arrays
    arrays = {}
    for var in var_list:
        arrays[var] = np.full(len(masses), np.nan)

    for i, mass in enumerate(masses):
        print(mass)

        # find last line
        with open(filenames[mass], "rb") as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
            last = f.readline()

        # get vars from last line
        for var, idx in var_list.items():
            arrays[var][i] = float(last.split()[idx])

    # add var columns to table
    for var in var_list:
        last_dats[var] = arrays[var]

    return last_dats


def get_max_data(mass, filenames, column):
    data = np.loadtxt(filenames[mass], usecols=(column,), unpack=True)
    return data[-1]


def get_expl_shok(masses, filenames):
    expl_dats = {}

    for mass in masses:
        expl_dats[mass] = {}
        time, ener, shok = np.loadtxt(filenames[mass], usecols=(0, 9, 11), unpack=True)
        index = np.max(np.where(ener < 1.e49))

        expl_dats[mass]['rmax'] = np.max(shok)
        expl_dats[mass]['texp'] = time[index]
        expl_dats[mass]['Eexp'] = np.max(ener)

    return expl_dats


def get_extra_ener(masses, filenames):
    extra_ener = {}
    for mass in masses:
        extra_ener[mass] = {}

        time, ener = np.loadtxt(filenames[mass], usecols=(0, 9), unpack=True)
        ener /= 1e51
        max_ener = ener[-1]

        extra_ener[mass]['finalEner'] = max_ener
        extra_ener[mass]['finalTime'] = time[-1]

        params = curve_fit(quadratic, time, ener)[0]

        extra_ener[mass]['asympTime'] = dquadzero(*params)
        extra_ener[mass]['asympEner'] = asymp_ener(*params)

    return extra_ener


def quadratic(x, a, b, c):
    return a + b * x + c * x ** 2


def dquadzero(b, c):
    return -b / (2 * c)


def asymp_ener(a, b, c):
    zero = dquadzero(b, c)
    return quadratic(zero, a, b, c)


# def compactMu(masses):
#    compMu = {}
#    for mass in masses:
#        compMu[mass] = {}
#        rad, mr, entr, ye = np.loadtxt('/mnt/research/SNAPhU/Progenitors2/s'
#               +str(mass)+'_presn.FLASH', usecols=(0,1,6,9), unpack=True, skiprows=11)
#        mr /= 2.e33 # put in solar masses
#        # First compute compactnesses
#        index = np.max(np.where(mr < 1.75))
#        compMu[mass]['xi1.75'] = 1.75/(rad[index]/1e8)
#        index = np.max(np.where(mr < 2.5))
#        compMu[mass]['xi2.5'] = 2.5/(rad[index]/1e8)
#        # On to T.Ertl parameters
#        #mu = np.diff(mr) / np.diff(rad) * 1e8
#        m1 = mr[np.where(entr >= 4.)][0]
#        r1 = rad[np.where(mr == m1)][0]
#        m2 = m1+0.3
#        r2 = rad[np.where(mr >= m2)][0]
#        mu = 0.3 / (r2-r1) * 1e8
#        #compMu[mass]['M4'] = np.interp(4., entr, mr)
#        #compMu[mass]['mu4'] = np.interp(4., entr[:-1], mu)
#        compMu[mass]['M4'] = m1
#        compMu[mass]['mu4'] = mu
#        # Now estimate He core mass
#        compMu[mass]['Mhe'] = np.max(mr[np.where(ye < 0.7)])
#    return compMu
