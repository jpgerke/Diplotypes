# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 13:00:02 2016

@author: justi
"""

import numpy as np
import pickle as pkl
from collections import Counter

############
# Calculates diplotype frequencies using views of array slices (no copying).
############

def stripit(h1, h2):
    res = []
    for x,y in zip(h1, h2):
        u = ''.join(map(str, x))
        v = ''.join(map(str, y))
        dip = sorted([u,v])
        res.append(''.join(map(str, dip)))
    return Counter(res)

def dipfreqs(sim):
    assert sim.shape[1] % 2 == 0
    m = sim.shape[1]//2
    gam1 = slice(0, m, None)
    gam2 = slice(m, None, None)
    hap1 = sim[:, gam1]
    hap2 = sim[:, gam2]
    counts = []
    for k in range(0,m-1):
        dip1 = hap1[:, k:(k+2)]
        dip2 = hap2[:, k:(k+2)]
        counts.append(stripit(dip1, dip2))
    return counts

if __name__ == '__main__':
    countdict = {}
    dsets = 'biparental fourway'.split()
    for x in range(1,5):
        g = str(x)
        for d in dsets:
            dat = np.loadtxt('../data/' + d + '_41_' + g + '_.txt', dtype=int)
            countdict['_'.join([d,g])] = dipfreqs(dat)
    pkl.dump(countdict, open("../data/countdict.pkl", 'wb'))