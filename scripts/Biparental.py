# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 11:12:10 2016

@author: justi
"""

from Diplotype import Diplotype, Population
import sympy
from sympy.plotting import plot
import itertools as it
import numpy as np
import matplotlib.pyplot as plt

r = sympy.Symbol('r')
A = Diplotype(('AA', 'AA'))
B = Diplotype(('BB', 'BB'))
F1 = Diplotype.cross(A,B)
F2 = F1.cycle()
F3 = F2.cycle()
F4 = F3.cycle()

for x in (F2,F3,F4):
    print(x.marginal().simplify())


AAAA = F4.unphase()['AAAA']
lamb = sympy.lambdify(r, AAAA, "numpy")
plt.plot(np.arange(0,0.5,0.05), lamb(np.arange(0,0.5,0.05)))

plot(AAAA, xlim=(0,0.5), ylim=(0,1))
plot(F4.unphase()['ABAB'], xlim=(0,0.5), ylim=(0,0.5))
plot(F2.unphase()['ABAB'], xlim=(0,0.5), ylim=(0,0.5))
