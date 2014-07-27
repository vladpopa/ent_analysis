#!/usr/bin/env python

import numpy
import cPickle

# pklfile = cPickle.load(open('pkl/core_stats.pkl', 'rb'))
pklfile = cPickle.load(open('pkl/condensed_stats.pkl', 'rb'))

k = len(pklfile.keys())
n = len(pklfile['z'])

stats = numpy.zeros(1, dtype = [(key, 'f8', n) for key in pklfile])

for key in pklfile:
   stats[key] = pklfile[key]

# numpy.save('npy/core_stats.npy', stats)
numpy.save('npy/condensed_stats.npy', stats)