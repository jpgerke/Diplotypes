# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 11:12:10 2016

@author: justi
"""

import sympy
import Diplotype as dp
import numpy as np


### Karl Broman's equations from Table 1 of Genetics 2012 ###
r,k = sympy.symbols('r k')

halfk = (1/2)**k #exponential decay due to fixation by homozygosity
poly = (1 - 2*r + 2*r**2)**(k-1) #polynomial present in all but decaying
absorb = ( (1-2*r)**k) / (1+2*r) #specific to absorbing states
decay = (1-2*r)**(k-1) # specific to decaying states

eqn1 = 1/(2*(1+2*r)) - (1/2)**(k+1)*(2 - poly + absorb) #AA|AA

eqn2 = r/(1+2*r) - (1/2)**(k+1)*(2 - poly - absorb) #AA|BB

eqn3 = halfk * (1 - poly) #AA|AB

eqn4 = halfk * (poly + decay) #AA|BB

eqn5 = halfk * (poly - decay) #AB|BA

eqs = (eqn1, eqn2, eqn3, eqn4, eqn5)
full = eqn1*2 + eqn2*2 + eqn3*4 + eqn4 + eqn5

#validate a few examples
tol = 1e-8
assert abs(full.evalf(subs={r:0.2,k:2})) - 1.0 < tol
assert abs(full.evalf(subs={r:0.5,k:3})) - 1.0 < tol
assert abs(full.evalf(subs={r:0.0,k:2})) - 1.0 < tol
assert abs(full.evalf(subs={r:0.2,k:4})) - 1.0 < tol
#####

### Map Broman's equations to the unphased genotypes
eqmap = {'AAAA': eqn1,  #AA
         'BBBB': eqn1,  #BB
         'AABB': eqn2,  #AB
         'BBAA': eqn2,  #BA
         'AAAB': eqn3,  #AH
         'ABAA': eqn3,  #HA
         'ABBB': eqn3,  #HB
         'BBAB': eqn3,  #BH
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

if __name__ == '__main__':
    ##### Prove the Diplotype object gives the same results as the equations
    #centimorgan range from 0 to 100 (in Morgans)
    cM = np.arange(0,1,0.01)
    #convert to recombination probabilities
    rvec = dp.recfun(cM)
    
    pops = (F2, F3, F4)
    for x, pop in enumerate(pops):
        gen = x + 2
        eqns = pop.tofunc()
        for key, value in eqns.items():
            mytest = eqmap[key].subs({k: gen})
            lam = sympy.lambdify(r, mytest, "numpy")
            diffs = np.abs(value(rvec) - lam(rvec))
            out = "Gen {} Type {}: Max difference {}".format(gen, key, np.max(diffs))
            print(out)
            assert np.max(diffs) < tol
    #####
