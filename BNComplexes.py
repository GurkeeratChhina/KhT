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

import itertools as itertools
import numpy as np
import math
from random import choice
from KhT import *
from Tangles import *
from Complex import *
from Drawing import *
from Cobordisms import *
from time import time
from fractions import Fraction

def ToExponent(exponent):
    return str(exponent).translate(str.maketrans("-0123456789.", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹·"))

class BNobj(object):
    """A BNobject is a pair [idempotent,q,h,delta(optional)], where idempotent is either 0 (b=solid dot) or 1 (c=hollow dot). 
    """
    # gr is a list of [qu, hom] of the bigradings, where qu is the quantum grading, hom is the homological grading
    __slots__ = 'idem','q','h','delta'
    
    def __init__(self,idempotent,q,h,delta="default"):
        self.idem = idempotent
        self.q = q
        self.h = h
        if delta=="default":
            self.delta = q/2-h
        else:
            self.delta = delta
    
    def idem2dot(self):
        if self.idem==0:
            return "●"#b (solid dot)
        else: 
            return "○"#c (hollow dot)
    
    def BNobj2String(self,switch="idem",index=-1):        
        
        if "idem" in switch:
            idem=self.idem2dot()
        else:
            idem = ""
        
        if (index == -1) or ("index" not in switch):
            index = ""
        else:
            index = str(index)
        if "q" in switch:
            q="q"+ToExponent(self.q)
        else:
            q=""
        if "h" in switch:
            h="h"+ToExponent(self.h)
        else:
            h=""
        if "delta" in switch:
            if 2*self.delta % 2 == 0: # delta is an integer
                delta="δ"+ToExponent(round(self.delta))
            else:
                #delta="δ"+ToExponent(round(2*self.delta))+"'²"
                delta="δ"+ToExponent(self.delta)
            
        else:
            delta=""
        
        grading=q+h+delta
        
        if (grading == "") or (index == ""):
            return index+grading+idem
        else:
            return index+":"+grading+idem
    
    def shift_q(self,shift): #shift q, keep h fixed; create new object
        return BNobj(self.idem,self.q+shift,self.h,self.delta+shift/2)
        
    def shift_h(self,shift): #shift h, keep q fixed; create new object
        return BNobj(self.idem,self.q,self.h+shift,self.delta-shift)

def coeff_simplify(num,field=2):
    if field > 1:
        return num % field
    else: # This probably needs to be fixed for field=0
        return num 

class BNmor_alt(object):# work in progress
    """An element of Bar-Natan's algebra is a list of pairs [power,coeff]
    'power' is an integer, which determines the exponent of D (if positive) and the exponent of S (if negative)
    'coeff' is some non-zero integer (= coefficient in the base ring/field) # Alternatively, a Fraction object
    """
    #__slots__ = 'pairs'
    
    def __init__(self,S,D,I):
        self.S = np.array(S) # list of coefficients
        self.D = np.array(D) # list of coefficients
        self.I = I #coefficient
    
    def simplify_BNmor(self,field=2):
        """simplify algebra elements by omitting superflous zeros."""
        self.S=[coeff_simplify(i) for i in self.S]
        self.D=[coeff_simplify(i) for i in self.D]
        self.I=coeff_simplify(self.I)
        while self.S[-1] ==0:
            del self.S[-1]
        while self.D[-1] ==0:
            del self.D[-1]
        return self
    
    def __add__(self, other):
        
        #newS = [a+b for a,b in zip_longest(self.S,other.S)]
        if len(self.S) < len(other.S):
            newS = other.S.copy()
            newS[:len(self.S)] += self.S
        else:
            newS = self.S.copy()
            newS[:len(other.S)] += other.S
        
        #newD = [a+b for a,b in zip_longest(self.D,other.D)]
        if len(self.D) < len(other.D):
            newD = other.D.copy()
            newD[:len(self.D)] += self.D
        else:
            newD = self.D.copy()
            newD[:len(other.D)] += other.D
        
        return BNmor(newS,newD,self.I+other.I).simplify_BNmor()

    def __mul__(self, other):
        newSmatrix = np.tensordot(self.S,other.S,axes=0)
        newS= [sum(np.diagonal(A[:, ::-1],len(other.S)-index)) for index in range(len(self.S)+len(other.S)-1)]
        
        newDmatrix = np.tensordot(self.D,other.D,axes=0)
        newD= [sum(np.diagonal(A[:, ::-1],len(other.D)-index)) for index in range(len(self.D)+len(other.D)-1)]
        
        return BNmor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0]).simplify_BNmor()
    
class BNmor(object):
    """An element of Bar-Natan's algebra is a list of pairs [power,coeff]
    'power' is an integer, which determines the exponent of D (if positive) and the exponent of S (if negative)
    'coeff' is some non-zero integer (= coefficient in the base ring/field) # Alternatively, a Fraction object
    """
    __slots__ = 'pairs','field'
    
    def __init__(self,pairs,field=2):
        self.pairs = pairs
        self.field = field
    
    def simplify_BNmor(self,field=2):
        """simplify algebra elements by adding all coeffients of the same power of D or S, omitting those with coefficient 0. This is very similar to simplify_decos"""
        if self.pairs == []:
            self.pairs=[]
        def droplast(l):
            return l[:-1]
        def add_coeffs(iterable):
            coeff=0
            for x in iterable:
                coeff+=x[-1]
            return coeff_simplify(coeff,field)
        self.pairs = [power+[add_coeffs(grouped)] for power,grouped in groupby(sorted(self.pairs),droplast)]
        self.pairs = [x for x in self.pairs if x[-1]!=0]
        return self
    
    def __add__(self, other):
        if self.pairs == []:
            return other
        if other.pairs == []:
            return self
        return BNmor(self.pairs+other.pairs,self.field).simplify_BNmor(self.field)

    def __mul__(self, other):
        if (self.pairs == []) or (other.pairs == []):
            return ZeroMor
        return BNmor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0],self.field).simplify_BNmor(self.field)
    
    def is_identity(self):
        if len(self.pairs)!=1:
            return False
        elif self.pairs[0][0]!=0:
            return False
        elif self.pairs[0][1] in [1,-1]:
            return True
        else:
            return False
    
    def contains_D(self):
        return all([pair[0]<=0 for pair in self.pairs])==False
    
    def contains_S(self):
        return all([pair[0]>=0 for pair in self.pairs])==False
    
    def negative(self): #create new morphism
        return BNmor([[pair[0],(-1)*pair[1]] for pair in self.pairs],self.field)
    
    def BNAlg2String(self):
        string=""
        for pair in self.pairs:
            coeff = pair[1]
            if (string != "") & (coeff > 0):# add plus sign if the next coefficient is positive, except for the first summand
                string += "+"
            if coeff < 0: # add minus sign in any case
                string += "-"
                coeff = abs(coeff)
            
            if coeff != 0:# omit any summands with coefficient 0
            
                exponent=abs(pair[0])
                if exponent==1: # omit exponent 1 from the notation
                    exponent = ""
                else:
                    exponent= ToExponent(exponent)
                
                if coeff==1: # omit coefficients 1 and -1 from the notation
                    coeff = ""
                else:
                    coeff = str(coeff) + "·"
                
                if pair[0] > 0:# powers of D
                    string += coeff + "D" + exponent
                if pair[0] < 0:
                    string += coeff + "S" + exponent
                if pair[0] == 0:
                    string += coeff + "id"
        return string
        
class BNComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of BNobj as labels on the vertices.
        - A matrix of BNmor as the adjacency matrix.
        - A field: 0 means Q, 1 means Z, prime p>1 means Z/pZ. 
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list """
    __slots__ = 'gens','diff','field'
    
    def __init__(self,gens,diff,field):
        self.gens = gens
        self.diff = np.array(diff)
        self.field = field
        self.toField(field)
    
    def toField(self,field):
        self.field = field
        self.diff = np.array([[mor.simplify_BNmor(field) for mor in row] for row in self.diff])
    
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
        for i,row in enumerate(squared):
            for j,mor in enumerate(row):
                mor.simplify_BNmor(self.field)
                if mor.pairs != []:
                    print("!!!!!!!!!!!!!!!!!!")
                    print("ERROR: Found non-zero term "+mor.BNAlg2String()+" in d² in row "+str(i)+" and column "+str(j)+".")
                    print("!!!!!!!!!!!!!!!!!!")
                    raise Exception('Differential does not square to 0')
        
    def isotopy(self,start,end,alg,switch="safe"):
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        
        if (switch=="safe") & ((self.diff)[start,end].pairs != []):
            raise Exception('This isotopy probably does not preserve the chain isomorphism type. There is an arrow going in the opposite direction of the isotopy.')
        self.diff[end,:]+=[alg.negative()*element for element in self.diff[start,:]] # subtract all precompositions with the differential (rows of diff)
        self.diff[:,start]+=[element*alg for element in self.diff[:,end]] # add all postcompositions with the differential (columns of diff)
        #self.ValidMorphism()
            
    
    def isotopy_via_vector_end(self,end,vector): # unused function: this is actually *much* slower
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        self.diff+=np.outer([x.negative() for x in vector],self.diff[end,:]) # subtract all precompositions with the differential (rows of diff)
        self.diff[:,end]+=np.dot(self.diff,vector) # add all postcompositions with the differential (columns of diff)
    
    def isolate_arrow(self,start, end, alg):
        """ Try to isotope away all arrows with the same start or the same end as 'arrow'.
        'arrow' is a list [start, end, alg], where 'alg' is a BNmor with a single entry.

        """
        def inverse(num): #this only works over a field;could replace this by a suitable function for positive characteristic
            if self.field == 0:
                return Fraction(1)/num
            elif self.field == 1:
                if num in [1, -1]:
                    return num
                else: 
                    raise Exception("Can't invert over Z")
            else: #taken from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
                s = 0
                S = 1
                r = self.field
                R = num
                while r != 0:
                    q = R // r
                    R, r = r, R - q * r
                    S, s = s, S - q * s
                #print((S*num)%self.field) # should be 1 if computed correctly
                return S
                
        face=alg.pairs[0][0]
        inverse_coeff=inverse(alg.pairs[0][1])
        
        def find_isotopy(index,entry):# Note: We are assuming here that there is at most one label to remove. 
            for pair in entry.pairs:
                if (pair[0]*face>0):
                    if ((self.diff)[index,end].pairs == []) & (index != end):
                        return BNmor([[pair[0]-face,pair[1]*inverse_coeff]],self.field)
            return ZeroMor
        
        # first remove all other arrows with the same start
        for index in range(len(self.gens)):
            if ((self.diff)[index,end].pairs == []) & (index != end): # check that isotopy is valid.
                for pair in self.diff[index,start].pairs:
                    if pair[0]*face>0:# same face
                        self.isotopy(end,index,BNmor([[pair[0]-face,pair[1]*inverse_coeff]],self.field),"unsafe")
        # attempt to speed up , but actually MUCH slower...
        #isotopy_vector=np.array([find_isotopy(index,entry) for index,entry in enumerate(self.diff[:,start])])
        #self.isotopy_via_vector_end(end,isotopy_vector)
        
        # secondly, remove all remaining arrows with the same end
        for index in range(len(self.gens)):
            if ((self.diff)[start,index].pairs == []) & (index != start): # check that isotopy is valid.
                for pair in self.diff[end,index].pairs:
                    if pair[0]*face>0:# same face
                        self.isotopy(index,start,BNmor([[pair[0]-face,(-1)*pair[1]*inverse_coeff]],self.field),"unsafe")
    
    def clean_up_once(self,SD):
        """ Simplify complex wrt the face D (1) or S (-1).
        """
        remaining=list(range(len(self.gens))) # list of unsimplified generators
        
        while remaining !=[]:
            #print(remaining)
            start_current=choice(remaining)
            #print(start_current)
            end_current = -1
            index_current = -1
            power_current = SD*math.inf
            
            changed=True
            
            while changed == True:
                
                changed=False
                for end in remaining: # find shortest arrow of the given face starting at start_current
                    for index, pair in enumerate(self.diff[end,start_current].pairs):
                        if pair[0]*SD > 0: # same face
                            if (power_current-pair[0])*SD > 0:
                                end_current = end
                                index_current = index
                                power_current = pair[0]
                                changed = True
            
                if end_current==-1: # if no element: subtract from remaining
                    end_current=start_current
                    start_current=-1
                
                changed=False
                for start in remaining: # try to find shorter arrow of the given face ending at end_current
                    for index, pair in enumerate(self.diff[end_current,start].pairs):
                        if pair[0]*SD > 0: # same face
                            if (power_current-pair[0])*SD > 0:
                                start_current = start
                                index_current = index
                                power_current = pair[0]
                                changed = True
                
            if (end_current == -1) or (start_current == -1): # if no element: subtract from remaining
                if end_current == -1:
                    remaining.remove(start_current)
                else:
                    remaining.remove(end_current)
                
            else: # isolate this shortest arrow
                self.isolate_arrow(start_current, end_current, BNmor([self.diff[end_current,start_current].pairs[index_current]],self.field))
                #print("start:", start_current,"end:", end_current, "isotopy:",self.diff[end_current,start_current].pairs[index_current])
                remaining.remove(start_current)
                remaining.remove(end_current)
    
    def is_looptype(self):# todo: This could be made much more efficient!!
        """Return 'True' is the complex is loop-type.
        """
        matrix_D=np.array([[entry.contains_D() for entry in row] for row in self.diff])
        matrix_S=np.array([[entry.contains_S() for entry in row] for row in self.diff])
        size=len(matrix_D)
        return  all([list(matrix_D[:,index]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_D[index,:]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_S[:,index]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_S[index,:]).count(True)<2 for index in range(size)])
    
    def clean_up(self,max_iter=200):
        """ Simplify complex alternatingly wrt D and S faces and stop after at most max_iter iterations. The default is 200 iterations. 
        """
        
        iter=0
        time0=time()
        time1=time0
        while iter<max_iter:
            if iter % 20 == 0: # check every now and then during the iteration, whether the complex is already loop-type.
                if self.is_looptype():
                    print("Clean-Up: Finished after "+str(iter)+" iteration(s) in "+str(round(time1-time0,1))+" second(s).")
                    break
            iter+=1
            self.clean_up_once(-1)# faces S
            self.clean_up_once(1) # faces D
            time2=time()
            print("iteration: "+str(iter)+" ("+str(round(time2-time1,1))+" sec)", end='\n')# testing how to monitor a process
            time1=time2
            #time.sleep(1)
        else:
            print("Clean-up: Terminated as a result of reaching maximum iteration value of", max_iter)
            
    def negative(self): #create new morphism
        return BNmor([[pair[0],(-1)*pair[1]] for pair in self.pairs],self.field)
    
    def shift_h(self,shift): # create new complex
        new_gens = [gen.shift_h(shift) for gen in self.gens]
        if shift % 2 == 0:
            new_diff = self.diff
        elif shift % 2 == 1:
            new_diff = [[alg.negative() for alg in row] for row in self.diff]
        else:
            raise Exception('Why are you trying to shift homological grading by something other than an integer? I cannot do that!')
        return BNComplex(new_gens,new_diff)
    
    def shift_q(self,shift): # modify the complex
        self.gens = [gen.shift_q(shift) for gen in self.gens]
    
    def cone(self,Hpower):
        
        shifted_complex = self.shift_h(1)
        shifted_complex.shift_q(2*Hpower)
        
        zero_matrix=np.array([[ZeroMor for i in range(len(self.gens))] for i in range(len(self.gens))])
        
        def fill_diagonal(i,j):
            if i==j:
                return BNmor([[Hpower,1],[-2*Hpower,(-1)**Hpower]],self.field)
            else:
                return ZeroMor
        
        Hdiagonal=np.array([[fill_diagonal(i,j) for j in range(len(self.gens))] for i in range(len(self.gens))])
        
        new_diff = np.concatenate(\
                (np.concatenate((self.diff,zero_matrix),axis=1),\
                np.concatenate((Hdiagonal,shifted_complex.diff),axis=1)),axis=0)
        
        return BNComplex(self.gens+shifted_complex.gens,new_diff)

ZeroMor=BNmor([])# note that this morphism is has field=2 by default, but this information gets never used.
#ZeroMor=BNmor([],[],0)

def CobordismToBNAlg(cob,field=2):
    """ Convert a cobordism between two (1,3)-tangles into an element of BNAlgebra."""
    
    if cob.decos==[]:
        return ZeroMor
    
    if cob.front.top !=1 or cob.front.bot !=3:
        raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
    
    if len(cob.comps)==1:# saddle
        return BNmor([[-1-2*deco[0],(-1**deco[0])*deco[-1]] for deco in cob.decos if deco[1]==0],field).simplify_BNmor(field)
        
    if len(cob.comps)==2:# identity/dot cobordism
        i=find_first_index(cob.comps,contains_0)+1 #component with TEI 0
        j=3-i #component without TEI 0: i=1 => j=2, i=2 => j=1
        
        decos_reduced=[deco for deco in cob.decos if deco[i]==0]     # ignoore those decos with dots in component with TEI 0
        decos_no_dots=[deco for deco in decos_reduced if deco[j]==0] #decos without dots
        
        decos_DD=[[deco[0]+1,deco[-1]]   for deco in decos_reduced if deco[j]==1] # dot cobordism
        decos_id=[[0,deco[-1]]           for deco in decos_no_dots if deco[0]==0] # id contribution
        decos_DH=[[deco[0],deco[-1]]     for deco in decos_no_dots if deco[0]>0 ] # D contribution from H
        decos_SH=[[-2*deco[0],-deco[-1]] for deco in decos_no_dots if deco[0]>0 ] # SS contribution from H
        
        return BNmor(decos_DD+decos_id+decos_DH+decos_SH,field).simplify_BNmor(field)

def CLT2BNObj(clt):
    """Convert a (1,3)-tangle into one of the two idempotents of BNAlg."""
    if clt.top !=1 or clt.bot !=3:
        
        raise Exception("The cobordism to convert to an element of BNAlgebra is not between (1,3)-tangles.")
    elif clt.arcs[0]==3:
        return BNobj(0,clt.qgr,clt.pgr) #b
    elif clt.arcs[0]==1:
        return BNobj(1,clt.qgr,clt.pgr) #c

def CobComplex2BNComplex(complex,field=2):
    gens=[CLT2BNObj(clt) for clt in complex.elements]
    diff=[[CobordismToBNAlg(cob,field) for cob in row] for row in complex.morphisms]
    return BNComplex(gens,diff,field)

def BNObj2CLT(bnobj):
    arcs = []
    if bnobj.idem == 0:
        arcs = [3,2,1,0]
    else:
        arcs = [1,0,3,2]
    return CLT(1,3, arcs, bnobj.h, bnobj.q, bnobj.delta)

def BNAlg2Cob(morphism, sourceCLT, targetCLT):
    decos = []
    for pair in morphism.pairs:
        if pair[0] > 0 : 
            decos.append([pair[0]-1, 0, 1, pair[1]])
        elif pair[0] == 0: 
            decos.append([0, 0, 0, pair[1]])
        elif pair[0] == 0 %2: 
            power = Fraction(-1*pair[0], 2)
            decos.append([power, 0, 0, ((-1)** power) *pair[1]]) #(-H)^(n/2)
            decos.append([power - 1, 0, 1, ((-1)** (power -1)) *pair[1]]) #(-H)^(n/2-1) D
        else: 
            power = Fraction(-1*pair[0] -1, 2)
            decos.append([power, 0, ((-1)** power)*pair[1] ])#(-H)^((n-1)/2) S
    Cob = Cobordism(sourceCLT, targetCLT, decos)
    # print("COB", Cob.front.top, Cob.front.bot)
    Cob.ReduceDecorations()
    return Cob
    
def BNComplex2CobComplex(BNcomplex):
    elements = [BNObj2CLT(bnobj) for bnobj in BNcomplex.gens]
    morphisms = [[BNAlg2Cob(morphism, elements[source], elements[target]) \
                 for source, morphism in enumerate(row)] for target, row in enumerate(BNcomplex.diff)]
    return ChainComplex(elements, morphisms)
    

def DrawBNComplex(complex, filename,vertex_switch="index_qhdelta"):
    "draw a graph of for the BNcomplex"
    g = Graph()
    size = len(complex.gens)
    g.add_vertex(size)
    canvas_size = (800*math.sqrt(size), 600*math.sqrt(size))
    
    Vertex_colour = g.new_vertex_property("string") # "black" is the horizontal CLT (b) and "white" is the vertical CLT (c)
    Vertex_labelling = g.new_vertex_property("string")
    Position = g.new_vertex_property("vector<float>")
    for i, gen in enumerate(complex.gens):
        if gen.idem == 0:
            Vertex_colour[g.vertex(i)] = "black"
        else:
            Vertex_colour[g.vertex(i)] = "white"
        Vertex_labelling[g.vertex(i)]=gen.BNobj2String(vertex_switch,i)
    vprops =  {'text' : Vertex_labelling,\
               'color' : "black",\
               'fill_color' : Vertex_colour,\
               'font_size' : 20,\
               'size' : 20}
    
    Edge_labeling = g.new_edge_property("string")
    for i in range(size):
        for j in range(size):
            edge_mor=complex.diff[j][i]
            if edge_mor.pairs != []:
                g.add_edge(g.vertex(i), g.vertex(j))
                Edge_labeling[g.edge(i,j)] = edge_mor.BNAlg2String()
    eprops =   {'color' : "black",\
                'pen_width' : 4.0,\
                'text' : Edge_labeling,\
                'text_color' : "black",\
                'text_distance' : 10,\
                'font_weight' : cairo.FONT_WEIGHT_BOLD,\
                'marker_size' : 20,\
                'font_size' : 22}   
                
    #position = arf_layout(g, max_iter=0)
    position = sfdp_layout(g, max_iter=0)
    #Position[g.vertex(i)] = [50*(2*i+1), 200]
    
    graph_draw(g, pos=position, vprops=vprops, eprops=eprops, output_size=canvas_size, bg_color=[1,1,1,1],  output="Output/" + filename)

def PrettyPrintBNComplex(complex):
    """Print a complex in human readable form.
    """
    print("The generators:")
    print(pd.DataFrame({\
        " ": [gen.idem2dot() for gen in complex.gens],\
        "q": [gen.q for gen in complex.gens],\
        "h": [gen.h for gen in complex.gens],\
        "δ": [gen.delta for gen in complex.gens]
        },columns=[" ","q","h","δ"]))
    print("The differential:")
    print(tabulate(pd.DataFrame([[entry.BNAlg2String() for entry in row] for row in complex.diff]),range(len(complex.diff)),tablefmt="fancy_grid"))


# Claudius: I'll keep working on this list... 
#todo: implement recognition of local systems (optional)
#todo: implement Cancellation (optional)
#todo: add a way to convert BNComplexes into Complexes (optional; could be useful for twisting)
#todo: implement pairing theorem (just for fun!)

BNmor0 = BNmor([])


#####################
####### TESTS #######
#####################

def Test_TwoTwistTangle():
    b=CLT(1,3,[1,0,3,2],[2,3])
    c=CLT(1,3,[3,2,1,0],[1,0])

    complex1 = ChainComplex([b,c,c,b,c,c], \
                [[ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                 [Cobordism(c,b, [[0,0,1]]) ,ZeroCob ,ZeroCob ,ZeroCob, ZeroCob, ZeroCob],\
                 [ZeroCob ,Cobordism(c, c, [[0,0,1,1]]) ,ZeroCob ,ZeroCob, ZeroCob, ZeroCob],\
                 [Cobordism(b, b, [[1,0,0,1]]) ,ZeroCob ,ZeroCob ,ZeroCob, ZeroCob, ZeroCob],\
                 [ZeroCob,Cobordism(c, c, [[1,0,0,1]]) ,ZeroCob ,Cobordism(c, b, [[0,0,-1]]) , ZeroCob, ZeroCob],\
                 [ZeroCob ,ZeroCob ,Cobordism(c, c, [[1,0,0,1]]) ,ZeroCob, Cobordism(c, c, [[0,0,1,-1]]), ZeroCob]])
    BNComplex1 = CobComplex2BNComplex(complex1)
    DrawBNComplex(BNComplex1, "TwoTwistTangle_before_cleanup.svg","index_h")
    PrettyPrintBNComplex(BNComplex1)

    #BNComplex1.isolate_arrow(0,2,BNmor([[-7,-1]]))
    #BNComplex1.isolate_arrow(1,2,BNmor([[-5,-1]]))
    #BNComplex1.isotopy(1,2,BNmor([[0,1]]))
    #BNComplex1.clean_up_once(-1)
    BNComplex1.clean_up()
    DrawBNComplex(BNComplex1, "TwoTwistTangle_after_cleanup.svg","index_h")
    PrettyPrintBNComplex(BNComplex1)
    
def Test_SplittingCurve():
    b=CLT(1,3,[1,0,3,2],[2,3])
    c=CLT(1,3,[3,2,1,0],[1,0])
    
    complex1 = ChainComplex([b,c,c,b,c,c], \
                [[ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                 [Cobordism(c,b, [[0,0,1]]) ,ZeroCob ,Cobordism(c, c, [[0,0,1,1]]) ,ZeroCob, ZeroCob, ZeroCob],\
                 [ZeroCob, ZeroCob ,ZeroCob ,ZeroCob, ZeroCob, ZeroCob],\
                 [Cobordism(b, b, [[1,0,0,1]]) ,ZeroCob ,ZeroCob ,ZeroCob, ZeroCob, ZeroCob],\
                 [ZeroCob,Cobordism(c, c, [[1,0,0,1]]) ,ZeroCob ,Cobordism(c, b, [[0,0,-1]]) , ZeroCob, Cobordism(c, c, [[0,0,1,-1]])],\
                 [ZeroCob ,ZeroCob ,Cobordism(c, c, [[1,0,0,1]]) ,ZeroCob, ZeroCob, ZeroCob]])
    BNComplex1 = CobComplex2BNComplex(complex1)
    DrawBNComplex(BNComplex1, "SplittingCurve_before_cleanup.svg","index_h")
    PrettyPrintBNComplex(BNComplex1)

    #BNComplex1.isolate_arrow(0,2,BNmor([[-7,-1]]))
    #BNComplex1.isolate_arrow(1,2,BNmor([[-5,-1]]))
    #BNComplex1.isotopy(1,2,BNmor([[0,1]]))
    #BNComplex1.clean_up_once(-1)
    BNComplex1.clean_up()
    DrawBNComplex(BNComplex1, "SplittingCurve_after_cleanup.svg","index_h")
    PrettyPrintBNComplex(BNComplex1)
    
def Test_2m3pt():# (2,-3)-pretzel tangle
    
    BNComplex1 = BNComplex(\
        [BNobj(1,-12,-5), BNobj(1,-10,-4),\
         BNobj(0,-11,-4), BNobj(0,-9,-3), BNobj(0,-7,-2),\
         BNobj(0,-9,-3),  BNobj(0,-7,-2), BNobj(0,-5,-1), BNobj(1,-4,0)],\
         [[BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor([[1,1],[-2,1]]),BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor([[-1,1]]),BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor0,BNmor([[-1,-1]]),BNmor([[-2,1]]),BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor0,BNmor0,BNmor0,BNmor([[1,1]]),BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor0,BNmor0,BNmor([[1,1]]),BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0],\
          [BNmor0,BNmor0,BNmor0,BNmor([[1,-1]]),BNmor0,BNmor([[-2,1]]),BNmor0,BNmor0,BNmor0],\
          [BNmor0,BNmor0,BNmor0,BNmor0,BNmor([[1,1]]),BNmor0,BNmor([[1,1]]),BNmor0,BNmor0],\
          [BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor0,BNmor([[-1,1]]),BNmor0]])
    DrawBNComplex(BNComplex1, "2m3pt_redBN_before_cleanup.svg","index_qh")
    PrettyPrintBNComplex(BNComplex1)

    #BNComplex1.isolate_arrow(0,2,BNmor([[-7,-1]]))
    #BNComplex1.isolate_arrow(1,2,BNmor([[-5,-1]]))
    #BNComplex1.isotopy(1,2,BNmor([[0,1]]))
    #BNComplex1.clean_up_once(-1)
    BNComplex1.clean_up()
    DrawBNComplex(BNComplex1, "2m3pt_redBN_after_cleanup.svg","index_qh")
    PrettyPrintBNComplex(BNComplex1)
    # This is the arc invariant = reduced Bar-Natan homology of the (2,-3)-pretzel tangle
    
    BNComplex2 = BNComplex1.cone(1)
    PrettyPrintBNComplex(BNComplex2)
    DrawBNComplex(BNComplex2, "2m3pt_redKh_before_cleanup.svg","index_qh")
    
    BNComplex2.clean_up()
    DrawBNComplex(BNComplex2, "2m3pt_redKh_after_cleanup.svg","index_qh")
    # This is the figure-8 invariant = reduced Khovanov homology of the (2,-3)-pretzel tangle
    
    BNComplex3 = BNComplex1.cone(2)
    PrettyPrintBNComplex(BNComplex3)
    DrawBNComplex(BNComplex3, "2m3pt_Kh_before_cleanup.svg","index_qh")
    
    BNComplex3.clean_up()
    DrawBNComplex(BNComplex3, "2m3pt_Kh_after_cleanup.svg","index_qh")
    # This is the lovely invariant = unreduced Khovanov homology of the (2,-3)-pretzel tangle

#Test_2m3pt()
#Test_2m3pt()
#Test_TwoTwistTangle()
#Test_SplittingCurve()

