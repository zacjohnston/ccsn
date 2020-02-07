import os
import numpy as np
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


def read_last_lines(masses, filenames):
    last_dats = {}

    for mass in masses:
        last_dats[mass] = {}

        with open(filenames[mass], "rb") as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
            last = f.readline()

        last_dats[mass]['ener'] = float(last.split()[9])
        last_dats[mass]['shok'] = float(last.split()[11])
        last_dats[mass]['dens'] = float(last.split()[16])
        last_dats[mass]['Mpns'] = float(last.split()[20])
        last_dats[mass]['time'] = float(last.split()[0])

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
