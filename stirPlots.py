import numpy as np
import os
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.colors as co
from scipy.optimize import fsolve, curve_fit
from scipy.integrate import quad

# adapted from https://github.com/snaphu-msu/ecRateStudy


def getBounceTime(masses, filenames):
    bouncetimes = {}

    for mass in masses:
        bouncetimes[mass] = {}
        logname = filenames[mass][:-4]+".log"
    
        for line in open(logname, 'r'):
            if "Bounce!" in line:
                bouncetimes[mass]["tbounce"] = float(line.split()[1])
    
        print(mass, bouncetimes[mass]["tbounce"])
    
    return bouncetimes


def readLastLines(masses, filenames):
    lastDats = {}
    
    for mass in masses:
        lastDats[mass] = {}
    
        with open(filenames[mass], "rb") as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
            last = f.readline()
    
        lastDats[mass]['ener'] = float(last.split()[9])
        lastDats[mass]['shok'] = float(last.split()[11])
        lastDats[mass]['dens'] = float(last.split()[16])
        lastDats[mass]['Mpns'] = float(last.split()[20])
        lastDats[mass]['time'] = float(last.split()[0])
    
    return lastDats


def getMaxData(mass, filenames, column):
    data = np.loadtxt(filenames[mass],usecols=(column,),unpack=True)
    return data[-1]


def getExplShok(masses, filenames):
    explDats = {}
    
    for mass in masses:
        explDats[mass] = {}
        time, ener, shok = np.loadtxt(filenames[mass], usecols=(0,9,11), unpack=True)
        index = np.max(np.where(ener < 1.e49))
    
        explDats[mass]['rmax'] = np.max(shok)
        explDats[mass]['texp'] = time[index]
        explDats[mass]['Eexp'] = np.max(ener)
    
    return explDats


def getExtraEner(masses, filenames):
    extraEner = {}
    for mass in masses:
        extraEner[mass] = {}

        time, ener = np.loadtxt(filenames[mass], usecols=(0,9), unpack=True)
        ener /= 1e51
        maxEner = ener[-1]
        
        extraEner[mass]['finalEner'] = maxEner
        extraEner[mass]['finalTime'] = time[-1]
        
        params = curve_fit(quadratic, time, ener)[0]
        
        extraEner[mass]['asympTime'] = dquadzero(*params)
        extraEner[mass]['asympEner'] = asympEner(*params)

    return extraEner


# def compactMu(masses):
#    compMu = {}
#    for mass in masses:
#        compMu[mass] = {}
#        rad, mr, entr, ye = np.loadtxt('/mnt/research/SNAPhU/Progenitors2/s'+str(mass)+'_presn.FLASH', usecols=(0,1,6,9), unpack=True, skiprows=11)
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


def quadratic(x, a, b, c):
    return a + b*x + c*x**2


def dquadzero(a,b,c):
    return -b/2./c


def asympEner(a,b,c):
    return quadratic(dquadzero(a,b,c), a, b,c)