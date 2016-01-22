# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:34:14 2016

@author: justi
"""
from __future__ import print_function
from collections import abc
import itertools as it

class Diplotype(object):
    """Two locus diplotype"""
    
    def __init__(self, haplotypes):
        if not isinstance(haplotypes, abc.Sequence):
            raise TypeError("haplotypes for diplotype must be a sequence.")
        if len(haplotypes) != 2:
            raise ValueError("Two haplotypes must be given (diploid).")
        if not all(len(x)==2 for x in haplotypes):
            raise ValueError("Each haplotypes must have two loci.")
        self.origin = haplotypes
        self.maternal = haplotypes[0]
        self.paternal = haplotypes[1]
       
    #The Standard methods   
    def __repr__(self):
        class_name = type(self).__name__
        return '{}({})'.format(class_name, self.origin)
    
    def __str__(self):
        return "|".join([self.maternal, self.paternal])
        
    def __iter__(self):
        return (i for i in self.origin)
    
    def __eq__(self, other):
        if not issubclass(type(other), type(self)):
            return False
        else:
            return self.origin == other.origin
    
    def __getitem__(self, x):
        return self.origin[x]

    def __len__(self):
        return len(self.origin)
    
    def __contains__(self, x):
        return x in self.alleles
    
    #class specific attributes
    @property
    def alleles(self):
        return {x for y in self.origin for x in y}

    @property
    def locus1(self):
        return ''.join([x[0] for x in self.loci])
    
    @property
    def locus2(self):
        return ''.join([x[1] for x in self.loci])
    
    #here is where we actually do something
    def gametes(self):
        """
        Generate gametes from the diplotype
                
        Returns a tuple of possible gametes.
        Probabilities are:
        (1-r, 1-r, r, r)
        
        """
        
        #if only one origin is present it's easy
        if len(self.alleles)==1:
            return([self.maternal[0]] * 4)
        else:
            prods = [self.maternal,
                     self.paternal,
                     ''.join([self.maternal[0], self.paternal[1]]),
                     ''.join([self.maternal[1], self.paternal[0]])]
            return(prods)
            
    def selfmate(self):
        gams = self.gametes()
        return [[gams[x], gams[y]]for x,y in it.product(range(0,4),
                                                        repeat = 2)]
                
if __name__ == '__main__':        
    b = Diplotype(["AA", "BB"])
    c = Diplotype(["AA", "AA"])
    d = Diplotype(["AB", "CD"])
    d.selfmate()
