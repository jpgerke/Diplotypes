# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:34:14 2016

@author: justin gerke

Written in python 3.5, but designed to also work in python 2.7
"""
from __future__ import print_function
from collections import abc, defaultdict
import itertools as it
import numpy as np
import sympy
from sympy import Symbol

#reverse haldane
def M_to_r(M):
    return 0.5*(1 - np.exp(-2*M))
    
#reverse haldane function for numpy vectors
recfun = np.vectorize(M_to_r)

class Diplotype(object):
    """
        Phased Two-locus diplotype
    
        Contains the parental gametes and the probability of origin.
        Probabilities can be reassigned, but genotypes are read-only.
        
        Doctests:
        >>> d = Diplotype(("AB", "CD"))
        >>> d
        Diplotype(('AB', 'CD'))
        
        >>> str(d)
        'AB|CD'
        
        >>> d == Diplotype(("AB", "CD"))
        True
        
        >>> tuple(d.gametes())
        (('AB', -r/2 + 1/2), ('CD', -r/2 + 1/2), ('AD', r/2), ('BC', r/2))
        
        >>> d.selfmate()[0]
        Diplotype(('AB', 'AB'), (-r/2 + 1/2)**2)
        
        >>> r = Symbol("r")
        >>> e = Diplotype(("AB","CD"), probability=r**2)
        >>> e
        Diplotype(('AB', 'CD'), r**2)
        
        >>> tuple(e.gametes())[0]
        ('AB', r**2*(-r/2 + 1/2))      
        
        >>> tuple(e.gametes())[1]
        ('CD', r**2*(-r/2 + 1/2))
        
        >>> e == d
        False
        
        >>> e == Diplotype(('AB','CD'), probability=r**2)
        True
        
        >>> A = Diplotype(('AA', 'AA'))
        >>> B = Diplotype(('BB', 'BB'))
        >>> F1 = Diplotype.cross(A,B)
        >>> sympy.simplify(F1.unphase()['ABAB'])
        1
    """
    
    def __init__(self, haplotypes, probability=None):
        if not isinstance(haplotypes, abc.Sequence):
            raise TypeError("haplotypes for diplotype must be a sequence.")
        if len(haplotypes) != 2:
            raise ValueError("Two haplotypes must be given (diploid).")
        if not all(len(x)==2 for x in haplotypes):
            raise ValueError("Each haplotypes must have two loci.")
        if probability is not None:
            if not issubclass(type(probability), sympy.Expr):
                raise TypeError("probability must of a sympy expression")
            self.prob = probability
        else:
            self.prob = 1
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
        elif self.__origin == other.__origin and self.prob == other.prob:
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
    def alleles(self): #Haven't used this yet but keeping it around
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
        probs = ((1-r)/2, (1-r)/2, r/2, r/2)
        prods = (self.maternal,
                 self.paternal,
                 ''.join([self.maternal[0], self.paternal[1]]),
                 ''.join([self.paternal[0], self.maternal[1]]))
        return zip(prods, probs)
            
    def selfmate(self):
        """Return all possible combinations of gametes"""
        zygotes = []
        for g,h in it.product(self.gametes(), repeat = 2):
            gamete_prob = g[1]*h[1]
            prob = self.prob * gamete_prob
            zygotes.append(Diplotype((g[0],h[0]), probability=prob))
        return Population(zygotes)

    def simplify(self):
        """
        Simplify probability if possible
        
        Slow. Only use once cycling is done.        
        """
        simp = sympy.simplify(self.__prob)   
        return Diplotype(self.__origin, simp)

    @classmethod
    def cross(cls, ind1, ind2):
        """Cross two diplotype objects"""
        zygotes = []
        for g,h in it.product(ind1.gametes(), ind2.gametes()):
            zygotes.append(cls((g[0],h[0]), probability=g[1]*h[1]))
        return Population(zygotes)                                                    

class Population(tuple):
    """An immutable sequence container for diplotypes"""
    
    def __new__(cls, x):
        return super().__new__(cls, (x))
    
    def combine(self):
        """
        Collect identical diplotypes
        
        Returns a new object of same class        
        """
        cls = type(self)
        t = defaultdict(int)
        for y in self:
           t[y._Diplotype__origin] += y.prob
        new = []
        for key, value in t.items():
            new.append(Diplotype(key, value))
        return cls(new)
        
    def unphase(self):
        """
        Combine while disregarding phase
        
        Returns a defaultdict         
        """
        t = defaultdict(int)
        for y in self:
            loc1 = sorted(y.locus1)
            loc2 = sorted(y.locus2)
            key = ''.join(loc1) + ''.join(loc2)
            t[key] += y.prob
        return t            

    def cycle(self):
        """Turn over a selfing population to the next generation."""
        #a generator of selfed Diplotypes from the Population
        nursery = (x.selfmate().combine() for x in self)
        
        #flatten into an iterator of all Diplotypes regardless of origin
        flat = it.chain.from_iterable(nursery)
        
        #construct a population from the iterator and combine like diplos
        return Population(flat).combine()

    def marginal(self):
        """
        Calculate marginal probability.

        Primarily for testing / debugging at this point.        
        """
        marg = 0
        for x in self:
            marg += x.prob
        return marg
        
    def tofunc(self):
        """Lambdify the equations for a population"""
        r = sympy.Symbol('r')
        newdict = {}
        for key, value in self.unphase().items():
            newdict[key] = sympy.lambdify(r, value, "numpy")
        return newdict         
     
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
    print(sympy.factor(probs['AA|AA']))
    Diplotype.cross(b,d)
    print("\nUnphased Genotype Frequencies for an F2 generation:")
    for x in sorted( b.selfmate().unphase().items() ):
        print(x)
    
