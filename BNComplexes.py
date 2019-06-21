import itertools as itertools
import numpy as np
import math
from KhT import *
from Tangles import *
from Cobordisms import *

class BNobj(object):
    """A BNobject is a pair [idempotent,gr], where 
    idempotent is either 0 (b) or 1 (c) 
    gr is a pair [quantum grading, homological grading] 
    ."""
    # gr is a list of [qu, hom] of the bigradings, where qu is the quantum grading, hom is the homological grading
    __slots__ = 'idempotent','gr'
    
    def __init__(self,idempotent,gr):
        self.idempotent
        self.gr = gr

class BNmor(object):
    """An element of Bar-Natan's algebra is a list of pairs [power,coeff]
    'power' is an integer, which determines the exponent of D (if positive) and the exponent of S (if negative)
    'coeff' is some non-zero integer (= coefficient in the base ring/field)
    """
    __slots__ = 'pairs'
    
    def __init__(self,pairs):
        self.pairs = pairs
    
    def simplify_BNmor(self):
        """simplify algebra elements by adding all coeffients of the same power of D or S, omitting those with coefficient 0. This is very similar to simplify_decos"""
        if self.pairs == []:
            self.pairs=[]
        def droplast(l):
            return l[:-1]
        def add_coeffs(l):
            return sum([element[-1] for element in l])
        self.pairs=[x for x in \
            [powers+[add_coeffs(list(grouped))] \
            for powers,grouped in groupby(sorted(self.pairs),droplast)] \
            if x[-1]!=0]
        return self
    
    def __add__(self, other):
        return BNmor(self.pairs+other.pairs).simplify_BNmor()

    def __mul__(self, other):
        return BNmor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0])
    
    def is_identity(self):
        if len(self.pairs)!=1
            return False
        elif self.pairs[0][0]!=0:
            return False
        elif self.pairs[0][1] in [1,-1]:
            return True
        else:
            return False

class BNComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of BNobj as labels on the vertices
        - A matrix of BNmor as the adjacency matrix.
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list """
    def __init__(self,gens,diff):
        self.gens = gens
        self.diff = diff
    
    def ValidMorphism(self):
        # return True
        length = len(self.gens)
        if len(self.diff) != length:
            raise Exception('Differential does not have n rows (where n is the number of elements in chain complex)')
        for i in self.diff:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of elements in chain complex)')
        for i in range(length):
            if len(self.diff[i][i].pairs)!=0:
                raise Exception('Differential has self loops')
        
        squared = np.tensordot(self.diff,self.diff, axes=(-2,-1))
        for i in flatten(squared):
            if i.ReduceDecorations() != []:
                raise Exception('Differential does not square to 0')
                
#todo: add a way to convert Complexes into BNComplexes
#todo: refactor the nice diagrammatic output for BNComplexes from 'Drawing.py'
#todo: add a way to convert BNComplexes into Complexes (optional; could be useful for twisting)
#todo: implement Clean-Up Lemma
#todo: implement Clean-Up Lemma to simplify at a given generator
#todo: implement randomized Clean-Upping (alternatingly D and S)
#todo: implement recognition of being loop-type 
#todo: implement recognition of local systems (optional)
#todo: implement Cancellation (optional)
#todo: implement pairing theorem (just for fun!)

BNmor0 = BNmor([])

print((BNmor([[0,1]])+BNmor([[0,-1],[2,24]])).pairs)
print((BNmor([[0,1]])*BNmor([[2,24]])).pairs)
print((BNmor([[-1,1]])*BNmor([[2,24]])).pairs)

