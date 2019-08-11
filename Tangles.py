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

import math

class CLT(object):
    """A crossingless tangle (CLT) is a matching of n+m ends (m=top,n=bot=bottom). 
    arcs is a list whose ith entry gives the value of the involution at the ith tangle end."""
    __slots__ = 'top','bot', 'total', 'arcs', 'pgr', 'qgr', 'dgr', 'pairs'
    
    def __init__(self,top,bot,arcs,pgr,qgr,dgr):
        self.top = top
        self.bot = bot
        self.total = int((top+bot)/2)# total number of arcs
        self.arcs = arcs
        self.pgr = pgr
        self.qgr = qgr
        self.dgr = dgr
        self.pairs = [pair for pair in [[i,arcs[i]] for i in range(top+bot)] if pair[0] < pair[1]] # testing alternative format for CLTs
    
    def shift_qhd(self,q,h,delta):
        self.pgr += h
        self.qgr += q
        self.dgr += delta
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
        return CLT(self.top+other.top,self.bot+other.bot,add(self.arcs,other.arcs),self.pgr+other.pgr, self.qgr + other.qgr, self.dgr+other.dgr) #TODO: check grading computations
    
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
        return CLT(self.top,other.bot,mul(self.arcs,other.arcs),self.pgr+other.pgr, self.qgr + other.qgr, self.dgr+other.dgr) #TODO: check grading computations

    def __eq__(self, other): # redefining eq so that the tangles don't have to be literally the same when composing cobordisms, but only need to have the same values
        """Check if two tangles are the same up to gradings.
        """
        #if self.top == other.top and self.bot == other.bot and self.arcs == other.arcs:
        #    if self.pgr == other.pgr and self.qgr == other.qgr and self.dgr == other.dgr:
        #        return True
        #    else:
        #        #print("tangles have different gradings, the first tangle has p,q,d:", self.pgr, self.qgr, self.dgr, "while the second has p,q,d:", other.pgr, other.qgr, other.dgr)
        #        return True
        #else:
        #    return False
        # Note: The above also checks gradings. However, it produces lots of messages during cancellation, since the homological grading changes. 
        return self.top == other.top and self.bot == other.bot and self.arcs == other.arcs
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
