# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 11:12:10 2016

@author: justi
"""

from Diplotype import Diplotype, Population
import sympy
import itertools as it
import numpy as np

r = sympy.Symbol('r')
A = Diplotype(('AA', 'AA'))
B = Diplotype(('BB', 'BB'))
F1 = Diplotype.cross(A,B)
F2 = F1.cycle()
F3 = F2.cycle()
F4 = F3.cycle()

for x in (F2,F3,F4):
    print(x.marginal().simplify())
