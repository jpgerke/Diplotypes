# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 11:12:10 2016

@author: justi
"""

import Diplotype as dp
import sympy
from sympy.plotting import plot
import numpy as np
import matplotlib.pyplot as plt


### Karl Broman's equations from Table 1 of Genetics 2012 ###
r = sympy.Symbol('r')
k = sympy.Symbol('k')

halfk = (1/2)**k
poly = (1 - 2*r + 2*r**2)**(k-1)
absorb = ( (1-2*r)**k) / (1+2*r) #specific to absorbing states
decay = (1-2*r)**(k-1) # specific to dec

eqn1 = 1/(2*(1+2*r)) - (1/2)**(k+1)*(2 - poly + absorb)

eqn2 = r/(1+2*r) - (1/2)**(k+1)*(2 - poly - absorb)

eqn3 = halfk * (1 - poly)

eqn4 = halfk * (poly + decay)

eqn5 = halfk * (poly - decay)

eqs = (eqn1, eqn2, eqn3, eqn4, eqn5)
full = eqn1*2 + eqn2*2 + eqn3*4 + eqn4 + eqn5

tol = 1e-8
assert full.evalf(subs={r:0.2,k:2}) - 1.0 < tol
assert full.evalf(subs={r:0.5,k:3}) - 1.0 < tol
assert full.evalf(subs={r:0.0,k:2}) - 1.0 < tol
assert full.evalf(subs={r:0.2,k:4}) - 1.0 < tol
#####

### Map Broman's equations to the unphased genotypes
eqmap = {'AAAA': eqn1,  #AA
         'BBBB': eqn1,  #BB
         'AABB': eqn2,  #AB
         'BBAA': eqn2,  #BA
         'AAAB': eqn3,  #AH
         'ABAA': eqn3,  #HA
         'ABBB': eqn3,  #HB
         'BBBA': eqn3,  #BH
         'ABAB': eqn4 + eqn5} #HH     
#####
         
### Diplotype frequencies from the model
A = dp.Diplotype(('AA', 'AA'))
B = dp.Diplotype(('BB', 'BB'))
F1 = dp.Diplotype.cross(A,B)
F2 = F1.cycle()
F3 = F2.cycle()
F4 = F3.cycle()
#####

#centimorgan range from 0 to 100 (in Morgans)
cM = np.arange(0,1,0.01)
rvec = dp.recfun(cM)
F2_eqs = F2.combine().tofunc()

#AAAA = F4.unphase()['AAAA']
#lamb = sympy.lambdify(r, AAAA, "numpy")
#plt.plot(np.arange(0,0.5,0.05), lamb(np.arange(0,0.5,0.05)))
#
#plot(AAAA, xlim=(0,0.5), ylim=(0,1))
#plot(F4.unphase()['ABAB'], xlim=(0,0.5), ylim=(0,0.5))
#plot(F2.unphase()['ABAB'], xlim=(0,0.5), ylim=(0,0.5))
