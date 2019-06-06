#!/usr/bin/env python3

#TO DO: 
#1. Implement Chain complexes
#2. Implement Chain complex operations such as flatten
#3. Implement delooping
#4. Implement gaussian elimination
#5. Implement tangle to CLT
#6. Test all of the above features
#7. Implement divide and conquer

# useful general functions that have nothing to do with this project
def find_first_index(mylist,function):
    """compute the index of the first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if function(v))

def find_first(mylist,function):
    """compute the index of the first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(v for v in mylist if function(v))

def indexQ(mylist,val):
    """compute the index of the first element in the list mylist which is equal to val. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if v==val)

def flatten(mylist):
    """analogue to Mathematica Flatten"""
    return [j for i in mylist for j in i]

#
# dictionary
#
# CLT= crossingless tangle
# TEI= tangle end index
# DS = alternate way of representing cobordisms of 4 ended tangles, using powers of D and powers of S

# load libraries
from itertools import groupby
from itertools import product
import itertools as itertools
import numpy as np
from graph_tool.all import *
import math
import cairo
from IPython.display import IFrame

###########################
###########################
###########################

class CLT(object):
    """A crossingless tangle (CLT) is a matching of n+m ends (m=top,n=bot=bottom). 
    arcs is a list whose ith entry gives the value of the involution at the ith tangle end."""
    # gr is a list of [qu, hom] of the bigradings, where qu is the quantum grading, hom is the homological grading
    __slots__ = 'top','bot', 'total', 'arcs', 'gr', 'pairs'
    
    def __init__(self,top,bot,arcs,gr):
        self.top = top
        self.bot = bot
        self.total = int((top+bot)/2)# total number of arcs
        self.arcs = arcs
        self.gr = gr 
        self.pairs = [pair for pair in [[i,arcs[i]] for i in range(top+bot)] if pair[0] < pair[1]] # testing alternative format for CLTs
    
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
            
    def __add__(self, other):
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
        return CLT(self.top+other.top,self.bot+other.bot,add(self.arcs,other.arcs),self.gr+other.gr)
    
    def __mul__(self, other):
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
        return CLT(self.top,other.bot,mul(self.arcs,other.arcs),self.gr+other.gr)
    
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
    def __init__(self,clt1,clt2,decos):
        self.front = clt1
        self.back = clt2
        self.comp = components(clt1,clt2)
        self.decos = decos
    
    def __add__(self, other):
        if self.decos == []: # Adding the zero cobordism
            return other
        if other.decos == []: # Adding the zero cobordism
            return self
        if self.front!=other.front or self.back!=other.back:# incompatible cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not compatible, because they do not belong to the same morphism space.')
        return Cobordism(self.front,self.back,simplify_decos(self.decos+other.decos))
        

    def __mul__(self, other):
        
        def partition_dC(dC,arcs):
            """find the finest partition of dC that is preserved by arcs. The output is a list of list of indices.
            Note: If dC is ordered wrt the lowest TEI, then so will be the output."""
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
        if y!=other.front: # incomposable cobordisms
            raise Exception('The cobordisms {}'.format(self)+' and {}'.format(other)+' are not composable; the first ends on {}'.format(y)+' are not composable; the first ends on {}'.format(other.clt1))
        #components of the fully simplified cobordism from x to z, ordered according to their smallest TEI
        dC=components(x,z)
        comps1=components(x,y)
        comps2=components(y,z)
        
        #finding the boundary component of each c of C 
        C=partition_dC(dC,y.arcs)# as lists dc of indices of components
        #print("dC=",dC)
        #print("C=",C)
        
        def find_comp(comp):
            return comp in C
        
        comp1_indices=[[j for j,comp in enumerate(comps1) if comp[0] in dc] for dc in dC]
        comp2_indices=[[j for j,comp in enumerate(comps2) if comp[0] in dc] for dc in dC]
        #print("comp1_indices=",comp1_indices)
        #print("comp2_indices=",comp2_indices)
        
        # |C_i intersect c|
        def c_i_cap_c(c_i,c_flat):
            #return sum(map(lambda s: s in c_flat, c_i))
            return sum([1 for i in c_i if i[0] in c_flat])
        genus=[1-int(0.5*(c_i_cap_c(self.comp,flatten([dC[j] for j in c]))\
                         +c_i_cap_c(other.comp,flatten([dC[j] for j in c]))\
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
            #print(l)
            return [sum([Hpower]+[i[0] for i in l])]+flatten([i[1:-1] for i in l])+[coeff*np.prod([i[-1] for i in l])]
        
        #print("dC = ",dC)
        decos=[]
        for e1 in self.decos:
            for e2 in other.decos:
                partial_decos=[]
                for i,c in enumerate(C):
                    r=sum([e1[index+1] for index in comp1_indices[i]]+[e2[index+1] for index in comp2_indices[i]])# number of dots on c
                    g=genus[i]# genus of the closure of c
                    IdcI=len(c)# number of boundary components of c
                    #print(",g=",g,",IdcI=",IdcI,",e1=",e1,",e2=",e2,",r=",r,",comp1_indices[i]=",comp1_indices[i],",comp2_indices[i]=",comp2_indices[i])
                    partial_decos.append(decos_from_c(g,r,IdcI))   
                    #print([decos_from_c(g,r,IdcI)])
                # decos_from_pair (e1,e2)
                #print("partial_decos = ",partial_decos)
                decos+=[combine_decos(l,e1[0]+e2[0],e1[-1]*e2[-1]) for l in itertools.product(*partial_decos)]
        
        #print("decos = ",decos)
        # reorder the dots according to the order of the components
        def reorder_dots(order):
            return [0]+[i[1]+1 for i in sorted(list(zip(order,range(len(order)))))]+[len(order)+1]
        # final result
        #print("simplify_decos(decos)",simplify_decos(decos))
        #print("C",C)
        #print("flatten(C)",flatten(C))
        #print("reorder_dots(flatten(C))",reorder_dots(flatten(C)))
        # PREVIOUSLY: return Cobordism(x,z,simplify_decos(decos)[:,reorder_dots(flatten(C))])
        #print('reorder dots')
        #print(reorder_dots(flatten(C)))
        #print(simplify_decos(decos))
        Output = Cobordism(x,z,[[deco[index] for index in reorder_dots(flatten(C))] for deco in simplify_decos(decos)])
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
            return len(x)==len(self.dots) and x==self.comp
        else:
            return False
        
    #def deg_safe(self):
        # check all summands in this linear combination to make sure the element is homogeneous.
        # FIXME
    
    # self is a linear combination of cobordisms, with each cobordism being a different list in decos.
    # ReduceDecorations sets a cobordism, decoration, with a dot on the top 0-th tangle end 
    # (which is chosen to be the basepoint of the cobordism) to be the zero cobordism, by removing decoration from decos
    # This also returns the resulting decorations
    def ReduceDecorations(self):
        ReducedDecorations = [decoration for decoration in self.decos if decoration[1] == 0]
        self.decos = ReducedDecorations
        return ReducedDecorations
    
def simplify_decos(decos):# ToDo: rewrite this without numpy
    if decos == []:
        return []
    decos=np.array(sorted(decos))
    """simplify decos by adding all coeffients of the same decoration, omitting those with coefficient 0."""
    # compute unique Hdots (unique_Hdot[0]) and how often each appears (unique_Hdot[1])
    unique_Hdot=np.unique(decos[:,:-1],return_counts=True,axis=0)
    # add all coefficients; we are assuming here that decos is ordered, otherwise it will not work!
    new_coeffs=[[sum(i)] for i in np.split(decos[:,-1],np.cumsum(unique_Hdot[1])[:-1])]
    # compute new decos
    decos=np.append(unique_Hdot[0],new_coeffs,axis=1)
    # return only those decorations whose coefficient is non-zero
    return (decos[decos[:,-1]!=0,:]).tolist()
    #return [x for x in [add_coeffs(list(g)) for k, g in groupby(sorted(decos),key = lambda s: s[0])] if x[1]!=0]

CLTA = CLT(1,1, [1,0], [0,0])
ZeroCob = Cobordism(CLTA,CLTA,[])

class ChainComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of CLTS as labels on the vertices
        - A a matrix of cobordisms as the adjacency matrix.
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list """
    def __init__(self,listofclt,morphisms):
        self.elements = listofclt
        self.morphisms = morphisms
    
    def ValidMorphism(self):
        length = len(self.elements)
        if len(self.morphisms) != length:
            raise Exception('Differential does not have n rows (where n is the number of elements in chain complex)')
        for i in self.morphisms:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of elements in chain complex)')
        for i in range(length):
            if self.morphisms[i][i].decos != []:
                raise Exception('Differential has self loops')
        
        squared = flatten(np.tensordot(self.morphisms,self.morphisms, axes=(-1,-2)))
        for i in squared:
            if i.decos != []:
                raise Exception('Differential does not square to 0')

def CobordismToDS(Cob):
    """ DS is a linear combination of morphisms that are powers of D or powers of S, represented as a list of lists
        an element of DS will be refered to a ds, which is a 3 element list
        the first element is D if the morphism is a power of D, and S if it is a power of S
        the second element is the power of the corresponding morphism
        the third element is the coefficient of the morphism
        Requires Cob to be a cobordism between 4 ended tangles """
    
    if Cob.front.total !=2 or Cob.back.total !=2:
        raise Exception("Cobordism to convert to DS is not between 4-ended tangles")
    DS = [] 
    for elem in Cob.decos:
        if len(elem) == 3: # elem is H^k S
            for ds in DS:
                if ds[0] == "S" and ds[1] == 2*elem[0] +1: #check if saddle is already there
                    ds[2] += elem[2]
                    break
            else:
                DS.append(["S", 2*elem[0]+1, elem[2]]) #otherwise add saddle
        elif len(elem) == 4 and elem[2] == 1: # elem is H^k D
            for ds in DS:
                if ds[0] == "D" and ds[1] == elem[0] +1: #check if dot is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["D", elem[0]+1, elem[3]]) #otherwise add dot
        elif len(elem) == 4 and elem[0] != 0: # elem is H^k times id
            for ds in DS:
                if ds[0] == "S" and ds[1] == 2*elem[0]: #check if saddle part is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["S", 2*elem[0], ((-1)**elem[0])*elem[3]]) #otherwise add saddle part
            for ds in DS:
                if ds[0] == "D" and ds[1] == elem[0]: #check if dot part is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["D", elem[0], elem[3]]) #otherwise add dot part
        else: #elem is id, elem should be [0, 0, 0, n]
            for ds in DS:
                if ds[0] == "S" and ds[1] == 0: #check if id is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["S", 0, elem[3]]) #otherwise add id
    return DS

def AddCap(Complex, i):
    """ Adds a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot """
    NewElements = []
    for clt in Complex.elements:
        newarcs = clt.arcs.copy()
        for j, x in enumerate(newarcs):
            if x >= clt.top+i:
                newarcs[j] += 2 # shifts arc ends that appear after the cap
        newarcs.insert(clt.top +i, clt.top+i) # inserting the cap
        newarcs.insert(clt.top +i, clt.top+i+1)
        NewCLT = CLT(clt.top, clt.bot+2, newarcs, clt.gr) 
        NewElements.append(NewCLT)
    length = len(Complex.elements)
    NewMorphisms = []
    for j ,row in enumerate(Complex.morphisms):
        NewRow=[]
        for k, cob in enumerate(row):
            DecosCopy = cob.decos.copy()
            magic_index = components(NewElements[j], NewElements[k]).index([NewElements[j].top+i, NewElements[j].top+i+1]) #computes the index of the new component of the cobordism corresponding to the identity sheet of the cap
            NewCob = Cobordism(NewElements[j], NewElements[k], [NewDeco.insert(magic_index,0) for NewDeco in DecosCopy])
            NewRow.append(NewCob)
        NewMorphisms.append(NewRow)
    return ChainComplex(NewElements, NewMorphisms)

def AddCupToCLT(clt, i):
    """ Here 0 <= i <= clt.bot
        Returns a list of 2 clt if there is a closed component
        otherwise returns a list of 1 clt """
    newElements = []
    if clt.arcs[clt.top +i] == clt.top+i+1: # adding the cup makes a closed component
        newarcs = clt.arcs.copy()
        newarcs.remove(clt.top +i) #removing the closed component
        newarcs.remove(clt.top+i+1)
        for j, x in enumerate(newarcs):
            if x >= clt.top+i:
                newarcs[j] -= 2 #shifting arc ends that appear after cup
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, [clt.gr[0]+1, clt.gr[1]]))
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, [clt.gr[0]-1, clt.gr[1]]))
    else: # adding the cup doesnt make closed components
        newarcs = clt.arcs.copy()
        leftend = newarcs[clt.top + i]
        rightend = newarcs[clt.top+i+1]
        newarcs[leftend] = rightend # connecting the two arcs via the cup
        newarcs[rightend] = leftend
        del newarcs[clt.top +i] #deleting the unneeded tangles
        del newarcs[clt.top +i]
        for j, x in enumerate(newarcs):
            if x >= clt.top+i:
                newarcs[j] -= 2 #shifting arc ends that appear after cup
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, clt.gr))
    return newElements

def AddCup(Complex, i):
    """ Here 0 <= i <= tangle.bot """
    newElements = []
    for clt in Complex.elements:
        newElements.extend(AddCupToCLT(clt, i))
    length = len(Complex.elements)
    NewMorphisms = []
    for j in range(length): # j is the target clt, TODO rewrite without range, use enumerate instead
        newRow = []
        nextRow = []
        for k in range(length): # k is the source clt
            if Complex.elements[k].arcs[clt.top +i] == clt.top+i+1 and Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # source and target are both closed, add 2 cobordisms to newRow and nextRow
                if Complex.morphisms[j][k].decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = []
                    newDecos2 = []
                    newDecos3 = []
                    newDecos4 = []
                    magic_index = components(Complex.elements[k], Complex.elements[j]).index([Complex.elements[k].top+i, Complex.elements[k].top+i+1]) #computes index of closed component
                    for deco in Complex.morphisms[j][k].decos:
                        decocopy = deco.copy()
                        del decocopy[magic_index + 1] # removing closed component
                        if deco[magic_index + 1] == 0: # no dot on component already prior
                            newDecos1.append(decocopy)
                            newDecos4.append(decocopy)
                        else: # is dot on component already
                            newDecos3.append(decocopy)
                            Hdecocopy = decocopy.copy()
                            Hdecocopy[0] += 1
                            newDecos4.append(Hdecocopy)
                    simplify_decos(newDecos1)
                    simplify_decos(newDecos2)
                    simplify_decos(newDecos3)
                    simplify_decos(newDecos4)
                    Cobordism1 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[j], i)[0], newDecos1)
                    Cobordism2 = Cobordism(AddCupToCLT(Complex.elements[k], i)[1], AddCupToCLT(Complex.elements[j], i)[0], newDecos2)
                    Cobordism3 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[j], i)[1], newDecos3)
                    Cobordism4 = Cobordism(AddCupToCLT(Complex.elements[k], i)[1], AddCupToCLT(Complex.elements[j], i)[1], newDecos4)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
                    nextRow.append(Cobordism3)
                    nextRow.append(Cobordism4)
            elif Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # source is open but target is closed, add 1 cobordism to each
                if Complex.morphisms[j][k].decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = [] #TODO: Fill out decos
                    newDecos2 = [] #TODO: Fill out decos
                    simplify_decos(newDecos1)
                    simplify_decos(newDecos2)
                    Cobordism1 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[j], i)[0], newDecos1)
                    Cobordism2 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[j], i)[1], newDecos2)
                    newRow.append(Cobordism1)
                    nextRow.append(Cobordism2)
            elif Complex.elements[k].arcs[clt.top +i] == clt.top+i+1: # source is closed but target is open, add 2 cobordisms to newRow only
                if Complex.morphisms[j][k].decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = [] #TODO: Fill out decos
                    newDecos2 = [] #TODO: Fill out decos
                    simplify_decos(newDecos1)
                    simplify_decos(newDecos2)
                    Cobordism1 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[j], i)[0], newDecos1)
                    Cobordism2 = Cobordism(AddCupToCLT(Complex.elements[k], i)[1], AddCupToCLT(Complex.elements[j], i)[0], newDecos2)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
            else: # source and target are both open, add 1 cobordism to newRow only
                if Complex.morphisms[j][k].decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = [] #TODO: Fill out decos
                    simplify_decos(newDecos1)
                    Cobordism1 = Cobordism(AddCupToCLT(Complex.elements[k], i)[0], AddCupToCLT(Complex.elements[k], i)[1], newDecos1)
                    newRow.append(Cobordism1)
        NewMorphisms.append(newRow)
        if Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # target is closed
            NewMorphisms.append(nextRow)
    return ChainComplex(newElements, NewMorphisms)
       
# graphical output for a crossingless tangle

def draw_tangle_ends(posx,posy,clt,h,ctx):
    ctx.set_font_size(0.40)
    ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
    
    #dots with labels
    for i in range(0,clt.top):
        ctx.move_to(i+posx,posy)
        ctx.line_to(i+posx,posy)
        s = '{}'.format(i)
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(i - width/2+posx,posy-0.15)
        ctx.show_text(s)
    
    for i in range(0,clt.bot):
        ctx.move_to(i+posx,h+posy)
        ctx.line_to(i+posx,h+posy)
        s = '{}'.format(i+clt.top)
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(i - width/2+posx,posy+h+0.35)
        ctx.show_text(s)
    
    ctx.set_source_rgb(0,0,0)
    ctx.set_line_width(0.1)
    ctx.stroke()
    
def draw_arcs(posx,posy,clt,h,ctx,rgb):
    
    #arcs connecting dots
    for i in clt.arcs_top():
        ctx.move_to(posx+i[0],posy+0)
        ctx.curve_to(posx+i[0],posy+abs(i[1]-i[0])/2,posx+i[1],posy+abs(i[1]-i[0])/2,posx+i[1],posy+0)
    
    for i in clt.arcs_bot():
        ctx.move_to(posx+i[0]-clt.top,posy+h)
        ctx.curve_to(posx+i[0]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h)
        
    for i in clt.arcs_mix():
        ctx.move_to(posx+i[0],posy+0)
        ctx.curve_to(posx+i[0],posy+abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h)
    
    ctx.set_source_rgb(*rgb)
    ctx.set_line_width(0.06)
    ctx.stroke()
    
def draw_dot_on_arc(arc,clt,h,ctx,deco_index):
    ctx.set_line_width(0.3)
    
    if arc[0]<clt.top:
        ctx.move_to(2+arc[0],(h+1)*deco_index+0.1)
        ctx.line_to(2+arc[0],(h+1)*deco_index+0.1)
    else:
        ctx.move_to(2+arc[0]-clt.top,h+(h+1)*deco_index-0.1)
        ctx.line_to(2+arc[0]-clt.top,h+(h+1)*deco_index-0.1)
    
    # THE CODE IS NOT OPERATIONAL YET. SO FAR, WE DRAW COLOURED DOTS ON TANGLE ENDS, BUT IDEALLY DRAW DOTS ON ARCS. 
    #arcs connecting dots
    #if arc not in clt.arcs():
    #    raise Exception('Proper dot-decoration of cobordisms failed.')
    # 
    # if arc in clt.arcs_top():
    #    ctx.move_to(posx+i[0],posy+0)
    #    ctx.curve_to(posx+i[0],posy+abs(i[1]-i[0])/2,posx+i[1],posy+abs(i[1]-i[0])/2,posx+i[1],posy+0)
    #
    #if arc in clt.arcs_bot():
    #    ctx.move_to(posx+i[0]-clt.top,posy+h)
    #    ctx.curve_to(posx+i[0]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h)
    #    
    #if arc in clt.arcs_mix():
    #    ctx.move_to(posx+i[0],posy+0)
    #    ctx.curve_to(posx+i[0],posy+abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h)
    #
    ctx.set_source_rgb(255,0,0)
    ctx.stroke()

def drawclt(clt,name):
    scale = 100
    w = max(clt.top,clt.bot)
    h = max([1]+[abs(i[1]-i[0]) for i in clt.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt.top) for i in clt.arcs_mix()])
    
    surface = cairo.PDFSurface(name+'.pdf',w*scale,(h+1)*scale)
    ctx = cairo.Context(surface)
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    draw_tangle_ends(0,0,clt,h,ctx)
    draw_arcs(0,0,clt,h,ctx,(0,0,0))
    
    return IFrame(name+'.pdf', width='100%', height='300')

def drawcob(cob,name):
    clt1=cob.front
    clt2=cob.back
    scale = 100
    w = max(clt1.top,clt1.bot)+2
    h = max([1]+[abs(i[1]-i[0]) for i in clt1.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt1.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt1.top) for i in clt1.arcs_mix()]+\
                 [abs(i[1]-i[0]) for i in clt2.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt2.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt2.top) for i in clt2.arcs_mix()])
    
    surface = cairo.PDFSurface(name+'.pdf',w*scale,(h+1)*scale*len(cob.decos))
    ctx = cairo.Context(surface)
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)

    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    
    for deco_index,deco in enumerate(cob.decos):
        #arcs connecting dots in clt2 blue ("beta")
        draw_arcs(2,(h+1)*deco_index,clt2,h,ctx,(0,0,255))
        #arcs connecting dots in clt1 red ("alpha")
        draw_arcs(2,(h+1)*deco_index,clt1,h,ctx,(255,0,0))
        
        # draw dots on components
        for dot_index,dot in enumerate(deco[1:-1]):
            if dot in [0,1]:
                if dot==1:
                    draw_dot_on_arc(cob.comp[dot_index][0:2],clt1,h,ctx,deco_index)
            else:
                raise Exception('Some components are decorated by more than one dot.')
        
        # draw coefficient + H power
        ctx.set_font_size(0.40)
        ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
        ctx.set_source_rgb(0, 0, 0)
        draw_tangle_ends(2,0+(h+1)*deco_index,clt1,h,ctx)
        
        s = '{}.'.format(deco[-1])+'H^{}'.format(deco[0])
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(0, h/2+(h+1)*deco_index)
        ctx.show_text(s)
    
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(0.1)
        ctx.stroke()
    
    return IFrame(name+'.pdf', width='100%', height='300')

def DrawFourEndedChainComplex(complex, filename):
    # For now assume that all CLT are 2-2
    for CLT in complex.elements:
        if CLT.top != 2 or CLT.bot !=2:
            raise Exception("Not a four ended tangle")
    
    # If a CLT is a 2-2 tangle, then the horizontal tangle is [1,0,3,2] and the vertical is [2,3,0,1]
    g = Graph()
    size = len(complex.elements)
    g.add_vertex(size)
    # print("size = " + str(size))
    Vertex_labeling = g.new_vertex_property("string") # "white" is the horizontal CLT and "black" is the vertical CLT
    Position = g.new_vertex_property("vector<float>")
    for i, clt in enumerate(complex.elements):
        if clt.arcs[0] == 1:
            Vertex_labeling[g.vertex(i)] = "white"
        elif clt.arcs[0] == 2:
            Vertex_labeling[g.vertex(i)] = "black"
        else: 
            print(complex.elements[i].arcs)
            print(complex.elements[i].arcs[0])
            raise Exception("CLT to label vertex is not a horizontal or vertical CLT")
    
        #Position[g.vertex(i)] = [50*(2*i+1), 200]
    
    # TODO: omit coefficients of +- 1, leaving only the sign
    SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    Edge_labeling = g.new_edge_property("string") # construct edge labels with linear combinations of powers of S and D
    for i, j in itertools.product(range(0, size), range(0,size)): #rewrite with enumerate
        if complex.morphisms[i][j].ReduceDecorations() != []:
            g.add_edge(g.vertex(i), g.vertex(j))
            Edge_labeling[g.edge(i,j)] = ""
            for ds in CobordismToDS(complex.morphisms[i][j]):
                if Edge_labeling[g.edge(i,j)] != "" and ds[2] > 0:
                    Edge_labeling[g.edge(i,j)] += "+"
                if ds[2] != 0:
                    Edge_labeling[g.edge(i,j)] += str(ds[2]) + "·" + ds[0] + str(ds[1]).translate(SUP)

    graph_draw(g, vertex_color = "black", vertex_fill_color = Vertex_labeling, vertex_size = 20,\
                edge_color = "black", edge_pen_width = 4.0, edge_text = Edge_labeling, edge_text_color = "black",\
                edge_text_distance = 10, edge_font_weight = cairo.FONT_WEIGHT_BOLD, edge_font_size = 22, \
                output_size=(1200, 400), output=filename)

# End of drawing code
   
# create basic crossingless tangles
def cup_alt(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cup"""
    return CLT(n+1,n-1,\
               [n+1+j for j in range(i-1)]+\
               [i,i-1]+\
               [n+j for j in range(i,n)]+\
               [j for j in range(i-1)]+\
               [j+2 for j in range(i-1,n-1)]\
               ,0)

def parallel(n):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cap"""
    return CLT(n,n,[n+j for j in range(n)]+[j for j in range(n)],0)
def cup(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cap"""
    return parallel(i-1)+CLT(2,0,[1,0],0)+parallel(n-i)
def cap(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cup"""
    return parallel(i-1)+CLT(0,2,[1,0],0)+parallel(n-i)


b=CLT(2,2,[1,0,3,2],[0,0])
drawclt(b,"b")      

c=CLT(2,2,[2,3,0,1],[0,0])
drawclt(c,"c")

Sbc=Cobordism(b,c,[[0,0,1]])
drawcob(Sbc,"Sbc")

Scb=Cobordism(c,b,[[0,0,1]])
drawcob(Scb,"Scb")
MinusScb = Cobordism(c,b,[[0,0,-1]])

#composition of morphisms
#print((Sbc*Scb).decos)
#print((Sbc*Scb).comp)
drawcob((Sbc*Scb),"SSbb")

#composition of morphisms
#print((Scb*Sbc).decos)
#print((Scb*Sbc).comp)
drawcob((Scb*Sbc),"SScc")

# addition of cobordisms
#Sbc=Cobordism(b,c,[[0,0,1]])
testcobordism=Cobordism(b,c,[[2,1,2],[0,0,4]])
drawclt(cup(10,3)*cap(10,2)*cup(10,4)*cup(8,6)+cap(2,1),"test1")
drawclt(cap(5,5),"test2")
drawcob(testcobordism,"testcobordism")

T1=CLT(2,4,[2,5,0,4,3,1],[0,0])
T2=CLT(2,4,[2,3,0,1,5,4],[0,0])
T3=CLT(4,4,[4,7,3,2,0,6,5,1],[0,0])
T4=CLT(4,4,[6,7,3,2,5,4,0,1],[0,0])
T5=CLT(4,4,[3,2,1,0,7,6,5,4],[0,0])

#print(components(b,c))
cob1=Cobordism(T1,T2,[[4,1,0,1]])
cob2=Cobordism(T1,T2,[[4,1,0,-3],[2,0,1,1],[1,1,1,19]])
cob3=Cobordism(T2,T1,[[2,0,1,-2]])
cob4=cob1+cob2
cob5=cob1*cob3
cob6=cob3*cob1
cob5.decos

drawclt(T1,"CLT1")
drawclt(T2,"CLT2")
drawclt(T3,"CLT3")
drawclt(T4,"CLT4")
drawclt(T5,"CLT5")
drawcob(cob1,"cob1")
cob1.ReduceDecorations()
drawcob(cob1, "cob1Reduced")
drawcob(cob2,"cob2")
drawcob(cob3,"cob3")
drawcob(cob4,"cob4")
drawcob(cob5,"cob5")
drawcob(cob6,"cob6")

complex1 = ChainComplex([b,c], [[ZeroCob, ZeroCob], [Sbc ,ZeroCob]])
complex1.ValidMorphism()

CobRightDotMinusLeftDotVertical = Cobordism(c,c, [[0,0,1,1],[0,1,0,-1]])
drawcob(CobRightDotMinusLeftDotVertical, "temporary1")

complex2 = ChainComplex([b,c,c], [[ZeroCob, Sbc, ZeroCob],[ZeroCob, ZeroCob, CobRightDotMinusLeftDotVertical],[ZeroCob, ZeroCob, ZeroCob]])
complex2.ValidMorphism()
drawcob(Sbc*CobRightDotMinusLeftDotVertical, "temporary2")
RightDc = Cobordism(c,c,[[0,0,1,1]])
drawcob(RightDc, "RightDc")
drawcob(Sbc*RightDc, "DottedSaddle")

drawcob(ZeroCob, "ZeroCob")

DrawFourEndedChainComplex(complex2, "complex2.png")
DrawFourEndedChainComplex(complex1, "complex1.png")

RightDcMinusH = Cobordism(c,c,[[0,0,1,1], [1,0,0,-1]])

complex3 = ChainComplex([b,c,c,c], [[ZeroCob, Sbc, ZeroCob, ZeroCob],[ZeroCob,ZeroCob, RightDc, ZeroCob],[ZeroCob, ZeroCob, ZeroCob, RightDcMinusH],[ZeroCob, ZeroCob, ZeroCob, ZeroCob]])
DrawFourEndedChainComplex(complex3, "complex3.png")
complex3.ValidMorphism()

complex4 = ChainComplex([c,b,c,c], [[ZeroCob, MinusScb, RightDc, ZeroCob], [ZeroCob, ZeroCob, ZeroCob, Sbc], [ZeroCob,ZeroCob, ZeroCob, Cobordism(c,c,[[0,0,0,1]])], [ZeroCob, ZeroCob, ZeroCob, ZeroCob]])
DrawFourEndedChainComplex(complex4, "complex4.png")

complex5 = ChainComplex([c,b,b,c,b,b], [[ZeroCob, Scb, ZeroCob, Cobordism(c, c, [[1,0,0,1]]), ZeroCob, ZeroCob],\
                                        [ZeroCob, ZeroCob, Cobordism(b,b, [[0,0,1,1]]), ZeroCob, Cobordism(b,b,[[1, 0,0, -1]]), ZeroCob],\
                                        [ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, Cobordism(b,b,[[1, 0,0,1]])],\
                                        [ZeroCob, ZeroCob, ZeroCob, ZeroCob, Scb, ZeroCob],\
                                        [ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, Cobordism(b,b, [[0,0,1,1]])],\
                                        [ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob]])
DrawFourEndedChainComplex(complex5, "complex5.png")

#complex1cap0 = AddCap(complex1, 0)

newarcs1 = T1.arcs.copy()
for j, x in enumerate(newarcs1):
    if x >= T1.top:
        newarcs1[j] += 2
newarcs1.insert(T1.top, T1.top)
newarcs1.insert(T1.top, T1.top+1)
NewT1 = CLT(T1.top, T1.bot+2, newarcs1, T1.gr)
drawclt(NewT1, "NewT1")

newarcs2 = T2.arcs.copy()
for j, x in enumerate(newarcs2):
    if x >= T2.top:
        newarcs2[j] += 2
newarcs2.insert(T2.top, T2.top)
newarcs2.insert(T2.top, T2.top+1)
NewT2 = CLT(T2.top, T2.bot+2, newarcs2, T2.gr)
drawclt(NewT2, "NewT2")

DecosCopy = cob2.decos.copy()
magic_index = components(NewT1, NewT2).index([NewT1.top, NewT1.top+1])
for j, NewDeco in enumerate(DecosCopy):
    DecosCopy[j].insert(magic_index + 1, 0)
NewCob = Cobordism(NewT1, NewT2, DecosCopy)
drawcob(NewCob, "NewCob")



BasicComplex = ChainComplex([CLT(1,1, [1,0], [0,0])], [[ZeroCob]])
BasicCap = AddCap(BasicComplex, 1)

drawclt(BasicCap.elements[0], "basiccap")
Double = AddCup(BasicCap, 0)
drawclt(Double.elements[0], "double")

def PrintComplexMorphismMatrix(Complex):
    for i in Complex.morphisms:
        Row = []
        for j in i:
            Row.append(len(j.decos))
        print(Row)

PrintComplexMorphismMatrix(Double)

tempcob1 = Cobordism(CLT(1,3, [1, 0, 3, 2], [0,0]), CLT(1,3, [1, 0, 3, 2], [0,1]), [[0, 0, 1, 1]])
tempcomplex = ChainComplex([CLT(1,3, [1, 0, 3, 2], [0,0]), CLT(1,3, [1, 0, 3, 2], [0,1])], [[ZeroCob, ZeroCob], [tempcob1, ZeroCob]])
tempcomplexwithcup = AddCup(tempcomplex, 1)
PrintComplexMorphismMatrix(tempcomplexwithcup)
print("Decorations at 3, 0", tempcomplexwithcup.morphisms[3][0].decos)
print("Decorations at 3, 1", tempcomplexwithcup.morphisms[3][1].decos)
print("Grading of element at 0", tempcomplexwithcup.elements[0].gr)
print("Grading of element at 1", tempcomplexwithcup.elements[1].gr)
print("Grading of element at 2", tempcomplexwithcup.elements[2].gr)
print("Grading of element at 3", tempcomplexwithcup.elements[3].gr)

# print(components(T1, T2))

# todo:
# done) multiplication in cobordism category (medium)
# 2) stacking operation for objects and morphisms which lands in clts and elementary cobordisms (hard)
# 3) cancellation algorithm (easy)

# proof of concept for matrix multiplication for matrices with customized algebra addition and multiplication\n",

# class testalg(object):
	# def __init__(self,x):
		# self.x = x
	# def __mul__(self, other):
		# return testalg(2*self.x*other.x)
	# # def __rmul__(self, other):
		# # return testalg(2*self.x*other.x)
	# def __add__(self, other):
		# return testalg(2+self.x+other.x)

# a=testalg(1)
# b=testalg(2)
# c=testalg(3)
# d=testalg(4)

# A=[[a,b],[c,d]]
# print((np.tensordot(A,A, axes=(-1,-2)))[0,0].x)
	
# B=[[1,2],[3,4]]
# C=[[1,1],[0,0]]
# print(np.tensordot(C,B, axes=(-1,-2)))
print("File executed successfully")