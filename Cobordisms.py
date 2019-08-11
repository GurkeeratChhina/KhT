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
    - a list decos of (row) vectors v_i in a numpy array of decorations (decos) of the form [Hpower,dot_1,...,dot_n,coeff]:
        - the first entry Hpower is a non-negative integer which records the power of H
        - the entries dot_<i> specify the number of dots (0 or 1) on the ith component.  
        - the last entry coeff is some non-zero integer (= coefficient in the base ring/field)
    - components: a list of lists of TEIs, which be long to the same component; 
      the order of the entries dot_<i> for i=1,...,n corresponds to the order of components. 
      WARNING: No assumptions are made about each list, ie the TEIs may appear in *any* order!
      (We might want to change this at some point)
    """
    __slots__ = 'front','back','comps','decos'
    
    def __init__(self,clt1,clt2,decos,comps="default"):
        self.front = clt1
        self.back = clt2
        if comps=="default":
            self.comps = components(clt1,clt2)
        else:
            self.comps = comps
        self.decos = decos
    
    def print(self,switch="short"):
        if self.decos==[]:
            return ""
        else:
            if switch == "old long":
                return [self.comps,self.decos]
            if switch == "long":
                table=[["H:"]+[deco[0] for deco in self.decos]]+\
                      [[comp]+[deco[i+1] for deco in self.decos] \
                              for i,comp in enumerate(self.comps)]+\
                      [["coeff:"]+[deco[-1] for deco in self.decos]]
                tablealt=[["H:",0,1],[[1,2],0,1],["coeff:",3,4]]
                return tabulate(table,tablefmt="plain")
            else:
                return len(self.decos)
    
    def __add__(self, other):
        if self.decos == []: # Adding the zero cobordism
            return other
        if other.decos == []: # Adding the zero cobordism
            return self
        if self.front!=other.front or self.back!=other.back:# incompatible cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not compatible, because they do not belong to the same morphism space.')
        def reorder_decos(old_comps,new_comps,decos):
            """convert decos ordered wrt 'old_comps' into decos ordered wrt 'new_comps'
            """
            if old_comps==new_comps:
                return decos
            else:
                new_order=[indexMemberQ(old_comps,comp[0]) for comp in new_comps]
                return [[deco[0]]+[deco[index+1] for index in new_order]+[deco[-1]] for deco in decos]
        return Cobordism(self.front,self.back,simplify_decos(self.decos+reorder_decos(other.comps,self.comps,other.decos)),self.comps)

    def __mul__(self, other):
        """The composition of two cobordism is computed by performing neck-cutting along push-offs of all boundary components of the composition and then removing all closed components.
        """
        
        def partition_new_comps(new_comps,arcs):
            """The finest partition of new_comps such that the endpoints of any arc in arcs belong to the same subset. The output is a list of list of indices."""
            partition=[]
            # remaining components in our iteration
            remaining=list(range(len(new_comps)))
            while len(remaining)>0:
                # pick the first of the remaining components as the nucleus of a new element of new_comps and remove from remaining
                nucleus = [remaining.pop(0)]
                # list of TEIs of nucleus
                nucleusL = new_comps[nucleus[0]]
                # find all TEIs of nucleus which arcs connects to a different component (which has to be in the set of remaining components)
                joins = [i for i in nucleusL if arcs[i] not in nucleusL]
                while len(joins)>0:
                    # for the first element of joins, find the component and remove it from the remaining components (strictly reducing len(remaining))
                    new=remaining.pop(find_first_index(remaining,lambda s: arcs[joins[0]] in new_comps[s]))
                    # add this component to the nucleus of the new element of new_comps 
                    nucleus.append(new) 
                    # record the TEIs of the component
                    nucleusL=nucleusL+new_comps[new]
                    # add the new TEIs to joins if arcs sends them outside the nucleus
                    joins=[i for i in joins if arcs[i] not in nucleusL]
                partition.append(nucleus)
            return partition
        
        if self.decos == [] or other.decos == []: #Multiplying by the zero cobordism
            return ZeroCob
        if self.back != other.front: # incomposable cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not composable; the first ends on',self.back.top, self.back.bot, self.back.arcs, self.back.pgr, self.back.qgr, 'and the second starts on', other.front.top, other.front.bot, other.front.arcs, other.front.pgr, other.front.qgr)
        
        comps1=self.comps
        comps2=other.comps
        
        # list of components (lists of TEIs) of the fully simplified cobordism from self.front to other.back (after neck-cutting), ordered according to their smallest TEI
        new_comps=components(self.front,other.back)
        # list of indices of components of new_comps which belong to the same component before neck-cutting
        old_comps_x=partition_new_comps(new_comps,self.back.arcs)
        # list of lists of TEIs that belong to the same component before neck-cuttting
        old_comps=[flatten([new_comps[index] for index in old_comp_x]) for old_comp_x in old_comps_x]
        # list of lists of indices of components in first/second cobordism that belong to the old_comp in old_comps
        comps1_x=[[j for j,comp in enumerate(comps1) if comp[0] in old_comp] for old_comp in old_comps]
        # list of lists of indices of components in first/second cobordism that belong to the old_comp in old_comps
        comps2_x=[[j for j,comp in enumerate(comps2) if comp[0] in old_comp] for old_comp in old_comps]
         
        def comp_genus(old_comp,old_comp_x):
            def intersection(comp):
                return sum([1 for i in comp if i[0] in old_comp])
            return 1-(intersection(comps1)+intersection(comps2)-len(old_comp)//2+len(old_comp_x))//2
        genus=[comp_genus(old_comp,old_comp_x) for old_comp,old_comp_x in zip(old_comps,old_comps_x)]# genus of the closure of old_comps
        
        def decos_from_old_comp(g,r,n):
            if r>0:# r>0 and v_c=1
                return [[g+r-1]+[1 for j in range(n)]+[1]]
            else:
                if g%2==0:# genus even
                    return [[g+n-sum(dots)-1]+list(dots)+[(-1)**(g+n-sum(dots)-1)]\
                            for dots in list(product([0,1], repeat=n))[:-1]]
                if g%2==1:# genus odd
                    return [[g+n-sum(dots)-1]+list(dots)+[((-1)**(g+n-sum(dots)-1))]\
                            for dots in list(product([0,1], repeat=n))[:-1]]\
                            +[[g-1]+[1 for i in range(n)]+[2]]
        
        def combine_decos(l,Hpower,coeff):
            return [sum([Hpower]+[i[0] for i in l])]+flatten([i[1:-1] for i in l])+[coeff*prod([i[-1] for i in l])]
        
        decos=[]
        for deco1 in self.decos:
            for deco2 in other.decos:
                partial_decos=[]
                for comp1_x,comp2_x,gen,old_comp_x in zip(comps1_x,comps2_x,genus,old_comps_x):
                    r=sum([deco1[index+1] for index in comp1_x]+[deco2[index+1] for index in comp2_x])# number of dots on old_comp
                    partial_decos.append(decos_from_old_comp(gen,r,len(old_comp_x))) 
                decos+=[combine_decos(l,deco1[0]+deco2[0],deco1[-1]*deco2[-1]) for l in product(*partial_decos)]
        
        Output = Cobordism(self.front,other.back,simplify_decos(decos),[new_comps[index] for index in flatten(old_comps_x)])
        Output.ReduceDecorations() # This kills any cobordism in the linear combination that has a dot on the same component as the basepoint
        return Output
    
    def deg(self): #no parameter self.dots
        degrees = [sum(deco[:-1]) for deco in self.decos]
        if len(set(degrees))==1:
            return len(self.decos[0])-2-self.front.total-2*degrees[0]
        elif len(set(degrees))==0:
            return 0
        else:
            raise Exception('The linear combination of cobordism is not homogeneous. The decos are '+str(self.decos)+'.')
    
    def homogeneousQ(self):
        return len(set([sum(deco[:-1]) for deco in self.decos]))<=1
    
    def check(self): #no parameter self.dots, and self.comps may have different ordering than generated by components()
        if self.front.top==self.back.top and self.front.bot==self.back.bot:
            x=components(self.front,self.back)
            return len(x)==len(self.dots) and x==self.comps
        else:
            return False
    
    def ReduceDecorations(self):
        """delete all decorations in a cobordism ('self') that have a dot in the component containing the basepoint (the TEI 0)."""
        ReducedDecorations = [deco for deco in self.decos if deco[find_first_index(self.comps,contains_0)+1] == 0]
        self.decos = ReducedDecorations
        return ReducedDecorations
    
    def isIsom(self):
        """checks if self is the identity cobordism or the negative identity cobordism"""
        # check if it is a single cobordism (and not a linear combination of several cobordisms or the zero cobordism)
        # checks if the number of components of the cobordism is exactly the number of arcs of the front CLT (and hence also the back CLT)
        # checks if there are no dots or powers of H (dots and Hpowers are non-negative)
        # checks if the coefficient is +-1. 
        # Note that if one of these conditions is not satified, then python does not evaluate the following ones and immediately returns False. 
        return (len(self.decos) == 1) and (len(self.comps) == self.front.total) and (sum(self.decos[0][:-1])==0) and (self.decos[0][-1] in [1, -1]) 

    def negative(self):
        """returns a cobordism that is exactly the same as self, but with all the coefficients in the linear combination multiplied by -1"""
        self.decos = [ deco[:-1] + [deco[-1]*-1] for deco in self.decos]
        return self
        #return Cobordism(self.front, self.back, newDecos, self.comps)
        
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
