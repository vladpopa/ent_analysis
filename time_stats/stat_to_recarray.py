#!/usr/bin/env python
import numpy as np
import cPickle

pklfile = cPickle.load(open('pkl/core_time_stats.pkl', 'rb'))
k = len(pklfile.keys())
n = len(pklfile['z'])
stats = np.zeros(1, dtype = [(key, 'f8', n) for key in pklfile])
for key in pklfile:
   stats[key] = pklfile[key]
np.save('npy/core_time_stats.npy', stats)

pklfile = cPickle.load(open('pkl/condensed_time_stats.pkl', 'rb'))
k = len(pklfile.keys())
n = len(pklfile['z'])
stats = np.zeros(1, dtype = [(key, 'f8', n) for key in pklfile])
for key in pklfile:
   stats[key] = pklfile[key]
np.save('npy/condensed_time_stats.npy', stats)

pklfile = cPickle.load(open('pkl/plume_time_stats.pkl', 'rb'))
k = len(pklfile.keys())
n = len(pklfile['z'])
stats = np.zeros(1, dtype = [(key, 'f8', n) for key in pklfile])
for key in pklfile:
   stats[key] = pklfile[key]
np.save('npy/plume_time_stats.npy', stats)