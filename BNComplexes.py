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
        self.idempotent = idempotent
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
        return BNmor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0]).simplify_BNmor()
    
    def is_identity(self):
        if len(self.pairs)!=1:
            return False
        elif self.pairs[0][0]!=0:
            return False
        elif self.pairs[0][1] in [1,-1]:
            return True
        else:
            return False
    
    def BNAlg2String(self):
        SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        string=""
        for pair in self.pairs:
            if (string != "") & (pair[1] > 0):# add plus sign if the next coefficient is positive, except for the first summand
                string += "+"
            if pair[1] != 0:# omit any summands with coefficient 0
                exponent=abs(pair[0])
                if exponent==1: # omit exponent 1 from notation
                    exponent = ""
                else:
                    exponent= str(exponent).translate(SUP)
                if pair[0] > 0:# powers of D
                    string += str(pair[1]) + "·" + "D" + exponent
                if pair[0] < 0:
                    string += str(pair[1]) + "·" + "S" + exponent
                if pair[0] == 0:
                    string += str(pair[1]) + "·" + "id"
        return string
        

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



def CobordismToBNAlg(cob):
    """ Convert a cobordism between two (1,3)-tangles into an element of BNAlgebra."""
    
    if cob.front.top !=1 or cob.front.bot !=3:
        raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
    
    if len(cob.comps)==1:# saddle
        return BNmor([[-1-2*deco[0],deco[-1]] for deco in cob.decos if deco[1]==0]).simplify_BNmor()
        
    if len(cob.comps)==2:# identity/dot cobordism
        i=find_first_index(cob.comps,contains_0)+1 #component with TEI 0
        j=3-i #component without TEI 0: i=1 => j=2, i=2 => j=1
        
        decos_reduced=[deco for deco in cob.decos if deco[i]==0]     # ignoore those decos with dots in component with TEI 0
        decos_no_dots=[deco for deco in decos_reduced if deco[j]==0] #decos without dots
        
        decos_DD=[[deco[0]+1,deco[-1]]   for deco in decos_reduced if deco[j]==1] # dot cobordism
        decos_id=[[0,deco[-1]]           for deco in decos_no_dots if deco[0]==0] # id contribution
        decos_DH=[[deco[0],deco[-1]]     for deco in decos_no_dots if deco[0]>0 ] # D contribution from H
        decos_SH=[[-2*deco[0],-deco[-1]] for deco in decos_no_dots if deco[0]>0 ] # SS contribution from H
        
        return BNmor(decos_DD+decos_id+decos_DH+decos_SH).simplify_BNmor()

def CLT2BNObj(clt):
    """Convert a (1,3)-tangle into one of the two idempotents of BNAlg."""
    if clt.top !=1 or clt.bot !=3:
        raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
    elif clt.arcs[0]==1:
        return BNobj(0,clt.gr) #b
    elif clt.arcs[0]==3:
        return BNobj(1,clt.gr) #c

def Complexes2BNComplexes(complex):
    gens=[CLT2BNObj(clt) for clt in complex.elements]
    diff=[[CobordismToBNAlg(cob) for cob in row] for row in complex.morphisms]
    return BNComplex(gens,diff)

# Claudius: I'll keep working on this list... 
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

b=CLT(1,3,[1,0,3,2],[2,3])
c=CLT(1,3,[3,2,1,0],[1,0])

alg=[[1,2],[-1,34],[0,9],[3,0],[-342,999]]
print(alg)

print(BNmor(alg).BNAlg2String())

print(CLT2BNObj(b).idempotent)
print(CLT2BNObj(c).idempotent)

print(CobordismToBNAlg(Cobordism(b,c,[[1,0,24]])).pairs)
print(CobordismToBNAlg(Cobordism(b,b,[[1,0,1,24]])).pairs)
print(CobordismToBNAlg(Cobordism(b,b,[[1,1,1,24]])).pairs)
print(CobordismToBNAlg(Cobordism(b,b,[[1,0,0,24]])).pairs)

#print((BNmor([[0,1]])+BNmor([[0,-1],[2,24]])).pairs)
#print((BNmor([[0,1]])*BNmor([[2,24]])).pairs)
#print((BNmor([[-1,1]])*BNmor([[2,24]])).pairs)

