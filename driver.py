import json
import os

# local
import stirPlots

# adapted from https://github.com/snaphu-msu/ecRateStudy

# masses = [9.0,9.25,9.5,9.75,
#           10.0,10.25,10.5,10.75,
#           11.0,11.25,11.5,11.75,
#           12.0,12.25,12.5,12.75,
#           13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9,
#           14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9,
#           15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9,
#           16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9,
#           17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9,
#           18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9,
#           19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9,
#           20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9,
#           21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9,
#           22.0, 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9,
#           23.0, 23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9,
#           24.0, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6, 24.7, 24.8, 24.9,
#           25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9,
#           26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9,
#           27.0, 27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9,
#           28.0, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9,
#           29.0, 29.1, 29.2, 29.3, 29.4, 29.5, 29.6, 29.7, 29.8, 29.9,
#           30.0,31,32,33,35,40,45,50,55,60,70,80,100,120]

models_path = '/Users/zac/projects/data/stir/run_sukhbold/run_14may19_a1.25'
masses = [12.0, 13.0, 19.0]

# --- Generate filenames ---
filenames = {}
for mass in masses:
    mpath = os.path.join(models_path, f'run_{mass:.1f}')
    fname = f'stir2_14may19_s{mass:.1f}_alpha1.25.dat'
    filenames[mass] = os.path.join(mpath, fname)

lastDats = stirPlots.readLastLines(masses, filenames)
with open('lastDats.json', 'w') as f:
    json.dump(lastDats, f)

# explDats = getExplShok(['aprox','lab'], mass, filenames)
# with open('explDatsEC.json', 'w') as f:
# 	json.dump(explDats, f)

# bounceTimes = getBounceTime(mass, filenames['lab'])
# with open('bounceTimesLab.json', 'w') as f:
# 	json.dump(bounceTimes, f)

# bounceTimes = getBounceTime(mass, filenames['aprox'])
# with open('bounceTimesAprox.json', 'w') as f:
#     json.dump(bounceTimes, f)

#extraEner = getExtraEner(alpha, mass, filenames)
#with open('extraEner.json', 'w') as f:
#    json.dump(extraEner, f)

#compMu = compactMu(mass)
#with open('compMu2.json', 'w') as f:
#    json.dump(compMu, f)

#pltLandscapeEner(alpha, mass, filenames)
