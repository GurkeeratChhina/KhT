# -*- coding: utf-8 -*-
# COPYRIGHT 2019 Gurkeerat Chhina, Claudius Zibrowius
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from itertools import product
from itertools import groupby
import math
from KhT import *
from Tangles import *

def clt_front_pairs(components):
    """the tangle at the front, in TEI pair notation, of a cobordism with components 'components'."""
    return [component[i:i+2] for component in components for i in range(0,len(component),2)]

def clt_back_pairs(components):
    """the tangle at the back, in TEI pair notation, of a cobordism with components 'components'."""
    return [(component+component[0:1])[i+1:i+3] for component in components for i in range(0,len(component),2)]

def clt_front_arcs(components):
    """the tangle at the front, in TEI arc notation, of a cobordism with components 'components'."""
    return arc_to_involution(clt_front_pairs(components))

def clt_back_arcs(components):
    """the tangle at the back, in TEI arc notation, of a cobordism with components 'components'."""
    return arc_to_involution(clt_back_pairs(components))

class Cobordism(object):
    """A cobordism from a crossingless tangle clt1=front to another clt2=back consists of 
    - clt1=front, (these may just be pointers?)
    - clt2=back, (these may just be pointers?)
    - components (alternative to clt1 and clt2)
    - a list decos of (row) vectors v_i in a numpy array of decorations (decos) of the form [Hpower,dot1,...,dotn,coeff]:
        - the first entry Hpower is a non-negative integer which records the pwer of H
        - the entries dot<i> specify the number of dots (0 or 1) on the ith component. 
            The components are ordered such that each component starts with the lowest index appearing in that component
            and all components are ordered by that first entry. 
        - the last entry coeff is some non-zero integer (= coefficient in the base ring/field)
    """
    def __init__(self,clt1,clt2,decos,comps="default"):
        self.front = clt1
        self.back = clt2
        if comps=="default":
            self.comps = components(clt1,clt2)
        else:
            self.comps = comps
        self.decos = decos
    
    def __add__(self, other):
        if self.decos == []: # Adding the zero cobordism
            return other
        if other.decos == []: # Adding the zero cobordism
            return self
        if self.front!=other.front or self.back!=other.back:# incompatible cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not compatible, because they do not belong to the same morphism space.')
        return Cobordism(self.front,self.back,simplify_decos(self.decos+other.decos),self.comps)


    def __mul__(self, other):
        
        def partition_dC(dC,arcs):
            """find the finest partition of dC that is preserved by arcs. The output is a list of list of indices."""
            partition=[]
            # remaining components in our iteration
            remaining=list(range(len(dC)))
            while len(remaining)>0:
                # pick the first of the remaining components as the nucleus of a new element of dC and remove from remaining
                nucleus = [remaining.pop(0)]
                # list of TEIs of nucleus
                nucleusL = dC[nucleus[0]]
                # find all TEIs of nucleus which arcs connects to a different component (which has to be in the set of remaining components)
                joins = [i for i in nucleusL if arcs[i] not in nucleusL]
                while len(joins)>0:
                    # for the first element of joins, find the component and remove it from the remaining components (strictly reducing len(remaining))
                    new=remaining.pop(find_first_index(remaining,lambda s: arcs[joins[0]] in dC[s]))
                    # add this component to the nucleus of the new element of dC 
                    nucleus.append(new) 
                    # record the TEIs of the component
                    nucleusL=nucleusL+dC[new]
                    # add the new TEIs to joins if arcs sends them outside the nucleus
                    joins=[i for i in joins if arcs[i] not in nucleusL]
                partition.append(nucleus)
            return partition
        
        if self.decos == [] or other.decos == []: #Multiplying by the zero cobordism
            return ZeroCob
        x=self.front
        y=self.back
        z=other.back
        if self.back != other.front: # incomposable cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not composable; the first ends on',self.back.top, self.back.bot, self.back.arcs, self.back.pgr, self.back.qgr, 'and the second starts on', other.front.top, other.front.bot, other.front.arcs, other.front.pgr, other.front.qgr)
        #components of the fully simplified cobordism from x to z, ordered according to their smallest TEI
        dC=components(x,z)
        comps1=self.comps
        comps2=other.comps

        #finding the boundary component of each c of C
        C=partition_dC(dC,y.arcs)# as lists dc of indices of components
        
        comp1_indices=[[j for j,comp in enumerate(comps1) if comp[0] in dc] for dc in dC]
        comp2_indices=[[j for j,comp in enumerate(comps2) if comp[0] in dc] for dc in dC]
        
        # |C_i intersect c|
        def c_i_cap_c(c_i,c_flat):
            #return sum(map(lambda s: s in c_flat, c_i))
            return sum([1 for i in c_i if i[0] in c_flat])
        genus=[1-int(0.5*(c_i_cap_c(self.comps,flatten([dC[j] for j in c]))\
                         +c_i_cap_c(other.comps,flatten([dC[j] for j in c]))\
                         -0.5*len(flatten([dC[j] for j in c]))\
                         +len(c))) for c in C]
        
        def decos_from_c(g,r,IdcI):
            if r>0:# r>0 and v_c=1
                return [[g+r-1]+[1 for j in range(IdcI)]+[1]]
            else:
                if g%2==0:# genus even
                    return [[g+IdcI-sum(dots)-1]+list(dots)+[(-1)**(g+IdcI-sum(dots)-1)]\
                            for dots in list(product([0,1], repeat=IdcI))[:-1]]
                if g%2==1:# genus odd
                    return [[g+IdcI-sum(dots)-1]+list(dots)+[((-1)**(g+IdcI-sum(dots)-1))]\
                            for dots in list(product([0,1], repeat=IdcI))[:-1]]\
                            +[[g-1]+[1 for i in range(IdcI)]+[2]]
        
        def combine_decos(l,Hpower,coeff):
            return [sum([Hpower]+[i[0] for i in l])]+flatten([i[1:-1] for i in l])+[coeff*prod([i[-1] for i in l])]
        
        decos=[]
        for e1 in self.decos:
            for e2 in other.decos:
                partial_decos=[]
                for i,c in enumerate(C):
                    r=sum([e1[index+1] for index in comp1_indices[i]]+[e2[index+1] for index in comp2_indices[i]])# number of dots on c
                    g=genus[i]# genus of the closure of c
                    IdcI=len(c)# number of boundary components of c
                    partial_decos.append(decos_from_c(g,r,IdcI))
                    
                decos+=[combine_decos(l,e1[0]+e2[0],e1[-1]*e2[-1]) for l in product(*partial_decos)]
        
        Output = Cobordism(x,z,simplify_decos(decos),[dC[index] for index in flatten(C)])
        Output.ReduceDecorations() # This kills any cobordism in the linear combination that has a dot on the same component as the basepoint
        return Output
    # def __rmul__(self, other):
        # #print('__rmul__')
        # return other
    
    def deg(self):
        return len(self.dots)-self.front.total-2*self.Hpower-2*sum(self.dots)
    
    def check(self):
        if self.front.top==self.back.top and self.front.bot==self.back.bot:
            x=components(self.front,self.back)
            return len(x)==len(self.dots) and x==self.comps
        else:
            return False
        
    #def deg_safe(self):
        # check all summands in this linear combination to make sure the element is homogeneous.
        # FIXME
    
    def ReduceDecorations(self):
        """delete all decorations in a cobordism ('self') that have a dot in the component containing the basepoint (the TEI 0)."""
        ReducedDecorations = [deco for deco in self.decos if deco[find_first_index(self.comps,contains_0)+1] == 0]
        self.decos = ReducedDecorations
        return ReducedDecorations
    
    def isIsom(self):
        if len(self.decos) != 1: 
            return False
        if len(self.decos[0]) != self.front.total +2:
            return False
        for x in self.decos[0][:-1]:
            if x != 0:
                return False
        
        if self.decos[0][-1] != 1 and self.decos[0][-1] != -1:
            return False
        return True

    def negative(self):
        newDecos = [ deco[:-1] + [deco[-1]*-1] for deco in self.decos]
        return Cobordism(self.front, self.back, newDecos, self.comps)
        
def simplify_decos(decos):
    """simplify decos by adding all coeffients of the same decoration, omitting those with coefficient 0."""
    if decos == []:
        return []    
    def droplast(l):
        return l[:-1]
    def add_coeffs(l):
        return sum([element[-1] for element in l])
    return [x for x in \
        [decos_without_coeff+[add_coeffs(list(grouped))] \
        for decos_without_coeff,grouped in groupby(sorted(decos),droplast)] \
        if x[-1]!=0]
        
CLTA = CLT(1,1, [1,0], 0, 0, 0)
ZeroCob = Cobordism(CLTA,CLTA,[])
