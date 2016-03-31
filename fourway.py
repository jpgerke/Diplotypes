# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 07:23:55 2016

@author: justi
"""

import sympy
import dill
import Diplotype as dp

r,k = sympy.symbols('r k')

#==============================================================================
# Cross Design
#
# Symetric double-cross hybrid followed by selfing.
#
# Parents:  AA, BB, CC, DD
#
# Intermediates:  AA|BB, CC|DD
#
# F1: AA|BB * CC|DD
#==============================================================================

intermediate_1 = dp.Diplotype.from_string("AA|BB")
intermediate_2 = dp.Diplotype.from_string("CC|DD")

F1 = dp.Diplotype.cross(intermediate_1, intermediate_2).combine()
F2 = F1.cycle().combine()
F3 = F2.cycle().combine()
F4 = F3.cycle().combine()
F4list = list(F4)
F4list.sort(key = lambda x: str(x))
geno_order = [str(x) for x in F4list]

#identify diplotypes with equal probabilities
eqsets = dict()
counter = list(range(0, len(F4list)))
while len(counter) > 0:
    mydip = F4list[counter[0]]
    myset = set([str(mydip)])
    for dip in F4list[1::]:
        if dip.prob == mydip.prob:
            myset |= set([str(dip)])
    eqsets[str(mydip)] = myset
    for x in reversed(counter):
        if str(F4list[x]) in myset:
            counter.remove(x)

for key,value in eqsets.items():
    print('{}: {}\n'.format(key, value))

#genotypes = {str(x) for x in F2}
#geno_order = list(sorted(genotypes))
#
#equivalents = [{'A', 'B'}, {'C', 'D'}]
#
#transiter = [(str.maketrans(''.join(sorted(x)),
#                            ''.join(sorted(x, reverse=True))))
#                            for x in equivalents]
#eqsets = dict()
#counter = list(range(0, len(geno_order)))
#while len(counter) > 0:
#    myset = set([geno_order[counter[0]]])
#    for trans in transiter:
#        myset |= {x.translate(trans) for x in myset}
#    eqsets[geno_order[counter[0]]] = myset
#    for x in reversed(counter):
#        if geno_order[x] in myset:
#            counter.remove(x)

#symmetries = [('AA|AA', 'CC|CC'),
#              ('AB|AB', 'CD|CD'),
#              ('AA|AC', 'CA|CC'),
#              ('AA|CA', 'AC|CC')]
#
#for x,y in symmetries:
#    eqsets[x] |= eqsets[y]
#    eqsets.pop(y, None)

states = eqsets
stateorder = list(states.keys())

#Construct a diagonal matrix as a correction factor
#for the number of diplotypes per collapsed state
corrections = [sympy.nsimplify(1/len(states[x])) for x in stateorder]
#>>> print(corrections)
#[1/2, 1/2, 1/4, 1, 1]
correction = sympy.zeros(len(corrections), len(corrections))
for x in range(0,len(corrections)):
    correction[x,x] = corrections[x]

incidence = sympy.zeros(len(geno_order), len(stateorder))
#Genotype to collapsed state incidence matrix
for i,g in enumerate(geno_order):
    for key, value in states.items():
        if g in value:
            j = stateorder.index(key)
            incidence[i,j] = 1
            continue

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

### Convert the dict to a matrix
trans_list = []
for x in stateorder:
    trans_list.append([transitions[x][y] for y in stateorder])
transmat = sympy.Matrix(trans_list)

P,D = transmat.diagonalize()
P = sympy.simplify(P)
Pin = sympy.simplify(P**-1)
Pin = P**-1
#raise to power k
for i in range(0,D.shape[0]):
    D[i,i] = D[i,i]**(k)

#create the vector of initial genotypes
F1dict = {str(x): x.prob for x in F1}
F1vec = [x.prob if str(x) in F1dict.keys() else 0 for x in F4list]
F1vec = sympy.simplify(sympy.Matrix(F1vec).T)
#compute equations
res = F1vec*incidence*P*D*Pin*correction

simple = sympy.simplify(res)
fourway_eqs = {key: value for key, value in zip(stateorder, simple)}
dill.dump(fourway_eqs, open("../data/fourway_eq.pkl", 'wb'))
dill.dump(states, open("../data/fourway_statedict.pkl", 'wb'))