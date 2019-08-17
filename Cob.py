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

from itertools import product, groupby
from tabulate import tabulate

import BNAlgebra
from aux import *

class obj(object):
    """A crossingless tangle (CLT) is a matching of n+m ends (m=top,n=bot=bottom). 
    arcs is a list whose ith entry gives the value of the involution at the ith tangle end."""
    __slots__ = 'top','bot', 'total', 'arcs', 'h', 'q', 'delta', 'pairs'
    
    def __init__(self,top,bot,arcs,h,q,delta):
        self.top = top
        self.bot = bot
        self.total = int((top+bot)/2)# total number of arcs
        self.arcs = arcs
        self.h = h
        self.q = q
        self.delta = delta
        self.pairs = [pair for pair in [[i,arcs[i]] for i in range(top+bot)] if pair[0] < pair[1]] # testing alternative format for CLTs
    
    def shift_qhd(self,q,h,delta):
        self.h += h
        self.q += q
        self.delta += delta
        return self
    
    def check(self):
        top = self.top
        bot = self.bot
        arcs = self.arcs
        
        # some test functions to make sure the input is valid
        def ValidMappingQ(arcs):
            return all(0<=element<len(arcs) for element in arcs)
        def FixPointFreeQ(arcs):
            return all(arcs[i]!=i for i in arcs)
        def InvolutionQ(arcs):
            return all(arcs[arcs[i]]==i for i in arcs)
        def Straighten(top,bot,arcs):
            """Convert an (n+m)-tangle into a (m+n,0)-tangle by rotating the bottom tangle ends to the right of the top ones."""
            def aux(index):
                if index<top:
                    return index
                else:
                    return 2*top+bot-index-1
            return [aux(element) for element in arcs[0:top]+arcs[top+bot:top-1:-1]]
        def CrossingFreeQ(arcs):
            s=Straighten(top,bot,arcs)
            return all(all(min(i,s[i])<j<max(i,s[i]) for j in s[min(i,s[i])+1:max(i,s[i])]) for i in s)
        
        # now apply the tests defined above and, if they succeed, define the arcs 
        if ValidMappingQ(arcs):
            if FixPointFreeQ(arcs) and InvolutionQ(arcs) and CrossingFreeQ(arcs):
                return True
            else:
                raise Exception('arcs {} is not, but should be a fix-point free involution of integers between 0 and n, where n is odd.\n'.format(arcs)\
                                 +'\n- fixpoint free? {}'.format(FixPointFreeQ(arcs))\
                                 +'\n- involution? {}'.format(InvolutionQ(arcs))\
                                 +'\n- crossingless? {}'.format(CrossingFreeQ(arcs)))
        else:
            raise Exception('arcs {} is not a list l of integers between 0 and length(l)-1, so it cannot be a well-defined fix-point free involution.'.format(arcs))
    
    def arcs_all(self):
        return [[i,self.arcs[i]] for i in range(self.top+self.bot) if i<self.arcs[i]]
    
    def arcs_top(self):
        return [[self.arcs[i],i] for i in range(self.top) if self.arcs[i]<min(i,self.top)]
    
    def arcs_bot(self):
        return [[i,self.arcs[i]] for i in range(self.top,self.top+self.bot) if self.arcs[i]>max(i,self.top)]
    
    def arcs_mix(self):
        return [[i,self.arcs[i]] for i in range(self.top) if self.arcs[i]>=self.top]
    
    def arcs_all(self):
        return [[i,self.arcs[i]] for i in range(self.top+self.bot) if i<self.arcs[i]]
            
    def __add__(self, other): #not currently used
        def add(m1,m2):
            def aux1(index):
                if index<self.top:
                    return index
                else:
                    return other.top+index
            def aux2(index):
                if index<other.top:
                    return self.top+index
                else:
                    return self.top+self.bot+index
            return [aux1(m1[i]) for i in range(self.top)]+[aux2(m2[i]) for i in range(other.top)]+[aux1(m1[self.top+i]) for i in range(self.bot)]+[aux2(m2[other.top+i]) for i in range(other.bot)]
        return obj(self.top+other.top,self.bot+other.bot,add(self.arcs,other.arcs),self.h+other.h, self.q + other.q, self.delta+other.delta) #TODO: check grading computations
    
    def __mul__(self, other): #not currently used
        def mul(m1,m2):
            # aux1 finds the end of a path through the tangle ending at the top of the first tangle
            def aux1(index):
                middle=m1[index]
                while middle>=self.top: # ie if middle does not lie on the top of the first tangle
                    middle=m2[middle-self.top]
                    if middle<other.top: # ie if middle does not lie on the bottom of the second tangle
                        middle=m1[middle+self.top]
                    else: # ie if middle lies on the bot of the second tangle
                        return middle+self.top-other.top
                # if middle lies on the top of the first tangle
                return middle  
            # aux2 finds the end of a path through the tangle starting at the bottom of the second tangle
            def aux2(index):
                middle=m2[index]
                while middle<other.top: # ie if middle does not lie on the bottom of the second tangle
                    middle=m1[middle+self.top]
                    if middle>=self.top: # ie if middle does not lie on the top of the first tangle
                        middle=m2[middle-self.top]
                    else: # ie if middle lies on the top of the first tangle
                        return middle
                # if middle lies on the top of the first tangle
                return middle+self.top-other.top 
            
            return [aux1(i) for i in range(self.top)]+[aux2(i+other.top) for i in range(other.bot)]
        return obj(self.top,other.bot,mul(self.arcs,other.arcs),self.h+other.h, self.q + other.q, self.delta+other.delta) #TODO: check grading computations

    def __eq__(self, other): # redefining eq so that the tangles don't have to be literally the same when composing cobordisms, but only need to have the same values
        """Check if two tangles are the same up to gradings.
        """
        #if self.top == other.top and self.bot == other.bot and self.arcs == other.arcs:
        #    if self.h == other.h and self.q == other.q and self.delta == other.delta:
        #        return True
        #    else:
        #        #print("tangles have different gradings, the first tangle has p,q,d:", self.h, self.q, self.delta, "while the second has p,q,d:", other.h, other.q, other.delta)
        #        return True
        #else:
        #    return False
        # Note: The above also checks gradings. However, it produces lots of messages during cancellation, since the homological grading changes. 
        return self.top == other.top and self.bot == other.bot and self.arcs == other.arcs

    def ToBNAlgebra(self):
        """Convert a (1,3)-tangle into one of the two idempotents of BNAlg."""
        if self.top !=1 or self.bot !=3:
            
            raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
        elif self.arcs[0]==3:
            return BNAlgebra.obj(0,self.q,self.h) #b
        elif self.arcs[0]==1:
            return BNAlgebra.obj(1,self.q,self.h) #c
        
class mor(object):
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
    
    def __repr__(self):
        clt1=self.front
        clt2=self.back
        return "mor[obj[{},{},{},{},{},{}],obj[{},{},{},{},{},{}],{},{}]".format(clt1.top,clt1.bot,clt1.arcs,clt1.h,clt1.q,clt1.delta,clt2.top,clt2.bot,clt2.arcs,clt2.h,clt2.q,clt2.delta,self.decos,self.comps)
    
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
        if other == 0: # Adding the zero cobordism
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
        decos=simplify_decos(self.decos+reorder_decos(other.comps,self.comps,other.decos))
        if decos == []:
            return 0
        return mor(self.front,self.back,decos,self.comps)
    
    def __radd__(self, other):
        if other == 0: # Adding the zero cobordism
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
        decos=simplify_decos(self.decos+reorder_decos(other.comps,self.comps,other.decos))
        if decos == []:
            return 0
        return mor(self.front,self.back,decos,self.comps)

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
        
        if other == 0: #Multiplying by the zero cobordism
            return 0
        if self.back != other.front: # incomposable cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not composable; the first ends on',self.back.top, self.back.bot, self.back.arcs, self.back.h, self.back.q, 'and the second starts on', other.front.top, other.front.bot, other.front.arcs, other.front.h, other.front.q)
        
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
        
        Output = mor(self.front,other.back,simplify_decos(decos),[new_comps[index] for index in flatten(old_comps_x)])
        
        return Output.ReduceDecorations()# This kills any cobordism in the linear combination that has a dot on the same component as the basepoint
    
    def __rmul__(self, other):
        if other == 0: #Multiplying by the zero cobordism
            return 0
        raise Exception('This case should not occur.')
    
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
        r=self
        ReducedDecorations = [deco for deco in self.decos if deco[find_first_index(self.comps,contains_0)+1] == 0]
        if ReducedDecorations == []:
            return 0
        else:
            self.decos = ReducedDecorations
            return self
    
    def isIsom(self):
        """checks if self is the identity cobordism or the negative identity cobordism"""
        # check if it is a single cobordism (and not a linear combination of several cobordisms or the zero cobordism)
        # checks if the number of components of the cobordism is exactly the number of arcs of the front CLT (and hence also the back CLT)
        # checks if there are no dots or powers of H (dots and Hpowers are non-negative)
        # checks if the coefficient is +-1. 
        # Note that if one of these conditions is not satified, then python does not evaluate the following ones and immediately returns False. 
        return (len(self.decos) == 1) and (len(self.comps) == self.front.total) and (sum(self.decos[0][:-1])==0) and (self.decos[0][-1] in [1, -1]) 

    def __neg__(self):
        """returns a cobordism that is exactly the same as self, but with all the coefficients in the linear combination multiplied by -1"""
        # mutate the objects:
        # self.decos = [ deco[:-1] + [deco[-1]*-1] for deco in self.decos]
        # return self
        # safer option, just as fast:
        newDecos = [ deco[:-1] + [deco[-1]*-1] for deco in self.decos]
        return mor(self.front, self.back, newDecos, self.comps)
    
    def ToBNAlgebra(self,field=2):
        """ Convert a cobordism between two (1,3)-tangles into an element of BNAlgebra."""
        
        if self.front.top !=1 or self.front.bot !=3:
            raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
        
        if len(self.comps)==1:# saddle
            return BNAlgebra.mor([[-1-2*deco[0],((-1)**deco[0])*deco[-1]] for deco in self.decos if deco[1]==0],field).simplify_mor(field)
            
        if len(self.comps)==2:# identity/dot cobordism
            i=find_first_index(self.comps,contains_0)+1 #component with TEI 0
            j=3-i #component without TEI 0: i=1 => j=2, i=2 => j=1
            
            decos_reduced=[deco for deco in self.decos if deco[i]==0]     # ignoore those decos with dots in component with TEI 0
            decos_no_dots=[deco for deco in decos_reduced if deco[j]==0] #decos without dots
            
            decos_DD=[[deco[0]+1,deco[-1]]   for deco in decos_reduced if deco[j]==1] # dot cobordism
            decos_id=[[0,deco[-1]]           for deco in decos_no_dots if deco[0]==0] # id contribution
            decos_DH=[[deco[0],deco[-1]]     for deco in decos_no_dots if deco[0]>0 ] # D contribution from H
            decos_SH=[[-2*deco[0],-deco[-1]] for deco in decos_no_dots if deco[0]>0 ] # SS contribution from H
            
            return BNAlgebra.mor(decos_DD+decos_id+decos_DH+decos_SH,field).simplify_mor(field)
        
def components(clt1,clt2):
    """components of an elementary cobordism between two tangles (assuming the same clt.top and clt.bot)."""
    allcomponents=[]
    done=[]
    for i in range(2*clt1.total):#find component for ith tangle end
        if i not in done:#continue if we don't already have a component; otherwise move on
            cur = i#current position
            component = []
            while cur not in done:
                done.append(cur)
                component.append(cur)
                cur=clt1.arcs[cur]
                done.append(cur)
                component.append(cur)
                cur=clt2.arcs[cur]
            allcomponents.append(component)
    return allcomponents

def arc_to_involution(pairs):
    """convert a list of pairs of TEIs into an involution."""
    return [x[1] for x in sorted(pairs+[pair[::-1] for pair in pairs])]

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

