# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 16:02:14 2016

@author: justi
"""
import sympy
import numpy as np
import Diplotype as dp
from Biparental import eqs


r,k = sympy.symbols('r k')

F2 = dp.Diplotype.from_string('AA|BB').selfmate().combine()
genotypes = sorted([str(x) for x in F2])
#>>> print(genotypes)
#['AA|AA', 'AA|AB', 'AA|BA', 'AA|BB',
# 'AB|AA', 'AB|AB', 'AB|BA', 'AB|BB',
# 'BA|AA', 'BA|AB', 'BA|BA', 'BA|BB',
# 'BB|AA', 'BB|AB', 'BB|BA', 'BB|BB']

#F1 Diplotype frequencies in order of 'genotypes'
F1vec = sympy.zeros(1, len(genotypes))
F1vec[0,genotypes.index('AA|BB')] = 1
#>>> print(F1vec)
#Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

#Define equivalent states
states = {
 'AA|AA': {'AA|AA', 'BB|BB'},
 'AB|AB': {'AB|AB', 'BA|BA'},
 'AA|AB': {'AA|AB', 'AB|AA', 'BA|BB', 'BB|BA', 'AA|BA', 'AB|BB', 'BA|AA', 'BB|AB'},
 'AA|BB': {'AA|BB', 'BB|AA'},
 'AB|BA': {'AB|BA', 'BA|AB'}}

#Set the desired order for the collapsed states
stateorder= ['AA|AA', 'AB|AB', 'AA|AB', 'AA|BB', 'AB|BA']

#Construct a diagonal matrix as a correction factor
#for the number of diplotypes per collapsed state
#In the state collapse, the dipotype order does not matter
corrections = []
for state in stateorder:
    dips = states[state]
    unordered = {tuple(sorted(x.split('|'))) for x in dips}
    corrections.append(sympy.nsimplify(1/len(unordered)))
#>>> print(correction)
#[1/2, 1/2, 1/4, 1, 1]
correction = sympy.zeros(len(corrections), len(corrections))
for x in range(0,len(corrections)):
    correction[x,x] = corrections[x]
#>>> print(correction)
#Matrix([[1/2, 0, 0, 0, 0],
#        [0, 1/2, 0, 0, 0],
#        [0, 0, 1/4, 0, 0],
#        [0, 0, 0,  1,  0],
#        [0, 0, 0, 0,   1]])


incidence = sympy.zeros(len(genotypes), len(stateorder))
#Genotype to collapsed state incidence matrix
for i,g in enumerate(genotypes):
    for key, value in states.items():
        if g in value:
            j = stateorder.index(key)
            incidence[i,j] = 1
            continue
# >>> print(incidence)
        
# 'AA|AA', 'AB|AB', 'AA|AB', 'AA|BB', 'AB|BA'    
# Matrix([[1, 0, 0, 0, 0],   'AA|AA'
#         [0, 0, 1, 0, 0],   'AA|AB'
#         [0, 0, 1, 0, 0],   'AA|BA'
#         [0, 0, 0, 1, 0],   'AA|BB'
#         [0, 0, 1, 0, 0],   'AB|AA'
#         [0, 1, 0, 0, 0],   'AB|AB'
#         [0, 0, 0, 0, 1],   'AB|BA'
#         [0, 0, 1, 0, 0],   'AB|BB'
#         [0, 0, 1, 0, 0],   'BA|AA'
#         [0, 0, 0, 0, 1],   'BA|AB'
#         [0, 1, 0, 0, 0],   'BA|BA'
#         [0, 0, 1, 0, 0],   'BA|BB'
#         [0, 0, 0, 1, 0],   'BB|AA'
#         [0, 0, 1, 0, 0],   'BB|AB'
#         [0, 0, 1, 0, 0],   'BB|BA'
#         [1, 0, 0, 0, 0]])  'BB|BB'

#Initialize the transition dictionary
transitions = {x: {y: 0 for y in states.keys()} for x in states.keys() }

####Populate the dictionary
for x in states.keys():
    #self a diplotype of state x
    pop = dp.Diplotype.from_string(x).selfmate().combine()
    #then for each possible state y
    for y in states.keys():
        #for each progeny in the pop
        for p in pop:
            #add the probability if the progeny is in the state
            if str(p) in states[y]:
                    transitions[x][y] = sympy.simplify(transitions[x][y] + p.prob)
#>>> for key, value in transitions.items():
#    print('{}:{}'.format(key, value))
#AB|AB:{'AB|AB': 1, 'AA|AB': 0, 'AB|BA': 0, 'AA|AA': 0, 'AA|BB': 0}
#AA|AB:{'AB|AB': 1/4, 'AA|AB': 1/2, 'AB|BA': 0, 'AA|AA': 1/4, 'AA|BB': 0}
#AB|BA:{'AB|AB': (r - 1)**2/2, 'AA|AB': 2*r*(-r + 1), 'AB|BA': (r - 1)**2/2,
#       'AA|AA': r**2/2, 'AA|BB': r**2/2}
#AA|AA:{'AB|AB': 0, 'AA|AB': 0, 'AB|BA': 0, 'AA|AA': 1, 'AA|BB': 0}
#AA|BB:{'AB|AB': r**2/2, 'AA|AB': 2*r*(-r + 1), 'AB|BA': r**2/2,
#       'AA|AA': (r - 1)**2/2, 'AA|BB': (r - 1)**2/2}
     
### Convert the dict to a matrix
trans_list = []
for x in stateorder:
    trans_list.append([transitions[x][y] for y in stateorder])
transmat = sympy.Matrix(trans_list)

#Diagonalization
P,D = transmat.diagonalize()
P = sympy.simplify(P)
Pin = sympy.simplify(P**-1)
#raise to power k
for i in range(0,D.shape[0]):
    D[i,i] = D[i,i]**(k-1)

#>>> print([D[i,i] for i in range(0,5)])
#[(1/2)**(k - 1), 1, 1, (-r + 1/2)**(k - 1), (r**2 - r + 1/2)**(k - 1)]

#A set of closed-form equations for a biparental cross is of the form:
# F1 frequencies * Diplotype incidence matrix *
# eigenvectors * eigenvalues**k-1 * eigenvectors**-1 * correction factor

#Result
#P_0*Z*P*D^k*P^-1*C
res = F1vec*incidence*P*D*Pin*correction

###Validation
#centimorgan range from 0 to 100 (in Morgans)
cM = np.arange(0,1,0.01)
#convert to recombination probabilities
rvec = dp.recfun(cM)

tol = 1e-8
for gen in range(2,5):
    for eq in range(0,5):
        test = res[0,eq].subs({k:gen})
        test = sympy.lambdify(r, test, "numpy")
        val = eqs[eq].subs({k:gen})
        val = sympy.lambdify(r, val, "numpy")
        diffs = np.max(np.abs(val(rvec) - test(rvec)))
        assert(diffs < tol)
       # print(diffs)
    
    
    
    