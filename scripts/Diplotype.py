# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:34:14 2016

@author: justin gerke

Written in python 3.5, but designed to also work in python 2.7
"""
from __future__ import print_function
from collections import abc, defaultdict
import itertools as it
import sympy as sym
from sympy import Symbol

class Diplotype(object):
    """
        Two locus diplotype
    
        Contains the parental gametes and the probability of origin.
        Immutable to make mapping safe. 
    """
    
    def __init__(self, haplotypes, probability=None):
        if not isinstance(haplotypes, abc.Sequence):
            raise TypeError("haplotypes for diplotype must be a sequence.")
        if len(haplotypes) != 2:
            raise ValueError("Two haplotypes must be given (diploid).")
        if not all(len(x)==2 for x in haplotypes):
            raise ValueError("Each haplotypes must have two loci.")
        if probability is not None:
            if not issubclass(type(probability), sym.Expr):
                raise TypeError("probability must of a sympy expression")
            self.__prob = probability
        else:
            self.__prob = 1
        self.__origin = haplotypes
        self.__maternal = haplotypes[0]
        self.__paternal = haplotypes[1]
       
    #The Standard methods   
    def __repr__(self):
        class_name = type(self).__name__
        if self.prob == 1:
            return '{}({})'.format(class_name, self.__origin)
        else:
            return '{}({}, {})'.format(class_name, self.__origin, self.prob)
    
    def __str__(self):
        return "|".join([self.maternal, self.paternal])
           
    def __eq__(self, other):
        if not issubclass(type(other), type(self)):
            return False
        elif self.__origin == other.__origin and self.__prob == other.__prob:
             return True
        else:
            return False
 
    def __contains__(self, x):
        return x in self.alleles
    
    #class specific attributes
    @property
    def maternal(self):
        """maternal getter"""
        return self.__maternal
    
    @property
    def paternal(self):
        """paternal getter"""
        return self.__paternal
        
    @property
    def prob(self):
        """prob getter"""
        return self.__prob
    
    @property
    def alleles(self):
        return {x for y in self.__origin for x in y}

    @property
    def locus1(self):
        return ''.join(x[0] for x in self.__origin)
    
    @property
    def locus2(self):
        return ''.join([x[1] for x in self.__origin])
    
    #here is where we actually do something
    def gametes(self):
        """
        Generate gametes from the diplotype
                
        Returns a tuple of possible gametes.
        Probabilities are:
        (1-r, 1-r, r, r)
        
        """    
        r = Symbol("r")    
        probs = (x*self.prob for x in ((1-r)/2, (1-r)/2, r/2, r/2))
        prods = (self.maternal,
                 self.paternal,
                 ''.join([self.maternal[0], self.paternal[1]]),
                 ''.join([self.maternal[1], self.paternal[0]]))
        return tuple(zip(prods, probs))
            
    def selfmate(self):
        """Return all possible combinations of gametes"""
        gams = self.gametes()
        zygotes = []
        for g,h in it.product(gams, repeat = 2):
            zygotes.append(Diplotype((g[0],h[0]), probability=g[1]*h[1]))
        return zygotes
                                                        
     
if __name__ == '__main__':        
    b = Diplotype(["AA", "BB"])
    c = Diplotype(["AA", "AA"])
    d = Diplotype(["AB", "CD"])
    x = d.selfmate()
    print("selfing products of AA|BB:")
    print('\n'.join([str(x) for x in x]))
    print('\n')
    probs = defaultdict(float)
    for y in b.selfmate():
        probs[str(y)] += y.prob
    print("grouped by genotype:")
    print(probs.items())
    print('\n')
    print("AA|AA probability factored:")
    print(sym.factor(probs['AA|AA']))
    
