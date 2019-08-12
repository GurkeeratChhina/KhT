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

from itertools import groupby
import numpy as np
import math
import cairo
from graph_tool.all import *
from random import choice
from time import time
from KhT import *
from fractions import Fraction
from BNAlgebra import *
from subprocess import run
from time import sleep


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
        def simplify(mor):
            if mor is 0:
                return 0
            return mor.simplify_BNmor(field)
        self.diff = np.array([[simplify(mor) for mor in row] for row in self.diff])
    
    def __repr__(self):
        return "BNComplex[{},{},{}]".format(self.gens,self.diff,self.field)
    
    def __str__(self):
        """string representation of a BNcomplex; 
        usage:
        >>> print(BNcomplex)
        """
        string=("The generators:\n")
        string+=str(pd.DataFrame({\
            " ": [gen.idem2dot() for gen in self.gens],\
            "q": [gen.q for gen in self.gens],\
            "h": [gen.h for gen in self.gens],\
            "δ": [gen.delta for gen in self.gens]
            },columns=[" ","q","h","δ"]))
        string+="\nThe differential:\n"
        string+=str(tabulate(pd.DataFrame([[entry.BNAlg2String() for entry in row] for row in self.diff]),range(len(self.diff)),tablefmt="fancy_grid"))
        return string
    
    def draw(self, filename,vertex_switch="index_qhdelta"):
        "draw a graph of for the BNcomplex"
        g = Graph()
        size = len(self.gens)
        g.add_vertex(size)
        
        scaling_factor=150 #this may need some adjusting
        if "index" in vertex_switch:
            scaling_factor+=100
        if "q" in vertex_switch:
            scaling_factor+=50
        if "h" in vertex_switch:
            scaling_factor+=50
        if "delta" in vertex_switch:
            scaling_factor+=50
            
        canvas_size = (scaling_factor*math.sqrt(size), scaling_factor*math.sqrt(size))# square canvas
        
        Vertex_colour = g.new_vertex_property("string") # "black" is the horizontal CLT (b) and "white" is the vertical CLT (c)
        Vertex_labelling = g.new_vertex_property("string")
        Position = g.new_vertex_property("vector<float>")
        for i, gen in enumerate(self.gens):
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
                edge_mor=self.diff[j][i]
                if edge_mor is not 0:
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
        
        position = arf_layout(g, max_iter=0)
        #position = sfdp_layout(g, max_iter=0)
        #Position[g.vertex(i)] = [50*(2*i+1), 200]
        
        graph_draw(g, pos=position, vprops=vprops, eprops=eprops, output_size=canvas_size, bg_color=[1,1,1,1],  output="Output/" + filename)
    
    def ValidMorphism(self):
        # return True
        length = len(self.gens)
        if len(self.diff) != length:
            raise Exception('Differential does not have n rows (where n is the number of gens in chain complex)')
        for i in self.diff:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of gens in chain complex)')
        for i in range(length):
            if self.diff[i][i] is not 0:
                raise Exception('Differential has self loops')
        
        squared = np.tensordot(self.diff,self.diff, axes=(-2,-1))
        for i,row in enumerate(squared):
            for j,mor in enumerate(row):
                if mor is not 0:
                    mor.simplify_BNmor(self.field)
                    if mor is not 0:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: Found non-zero term "+mor.BNAlg2String()+" in d² in row "+str(i)+" and column "+str(j)+".")
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Differential does not square to 0')
        
    def findIsom(self): 
        """Returns the location of the first isomorphism it finds
           If no isomorphism is found, returns None"""
        for targetindex, row in enumerate(self.diff):
            for sourceindex, morphism in enumerate(row):
                if (morphism is not 0) and (morphism.is_isomorphism()):
                    return [sourceindex, targetindex]
        return None
    
    def eliminateIsom(self, sourceindex, targetindex):
        """Mutates self by eliminating the isomorphism at the specified location, via the gaussian elimination lemma
           Note that this does not check that the cobodism specified is actually an isomorphism.
        """
        
        Max=max(targetindex,sourceindex)
        Min=min(targetindex,sourceindex)
        coeff = self.diff[targetindex,sourceindex] # coefficient of the cancelling arrow
        
        if (len(coeff.pairs)!=1) or (coeff.pairs[0][0]!=0):
            raise Exception('You cannot cancel this arrow!')
        
        
        del self.gens[Max] # eliminate source and target from list of generators
        del self.gens[Min] # eliminate source and target from list of generators
        
        out_source = np.delete(self.diff[:,sourceindex],[Min,Max],0) #arrows starting at the source, omiting indices targetindex and sourceindex
        in_target = np.delete(self.diff[targetindex],[Min,Max],0) #arrows ending at the target, omiting indices targetindex and sourceindex
        
        def neg(entry):
            if entry is 0:
                return 0
            return entry.negative(self.field,inverse(coeff.pairs[0][1],self.field))
        
        in_target=np.array([neg(entry) for entry in in_target])
        
        self.diff=np.delete(self.diff,[Min,Max],0) # eliminate rows of indices targetindex and sourceindex
        self.diff=np.delete(self.diff,[Min,Max],1) # eliminate columns of indices targetindex and sourceindex
        
        self.diff = self.diff+np.transpose(np.tensordot(in_target,out_source, axes=0)) # update differential

    def eliminateAll(self): #mutates self by eliminating isomorphisms as long as it can find one
        while True:
            index_to_eliminate = self.findIsom()
            if index_to_eliminate is None:
                break
            else:
                self.eliminateIsom(index_to_eliminate[0], index_to_eliminate[1])
       
    def isotopy(self,start,end,alg,switch="safe"):
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        
        if (switch=="safe") and ((self.diff)[start,end] is not 0):
            raise Exception('This isotopy probably does not preserve the chain isomorphism type. There is an arrow going in the opposite direction of the isotopy.')
        self.diff[end,:]+=[alg.negative(self.field)*element for element in self.diff[start,:]] # subtract all precompositions with the differential (rows of diff)
        self.diff[:,start]+=[element*alg for element in self.diff[:,end]] # add all postcompositions with the differential (columns of diff)
        #self.ValidMorphism()
            
    
    def isotopy_via_vector_end(self,end,vector): # unused function: this is actually *much* slower
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        def neg(x):
            if x is 0:
                return 0
            return x.negative(self.field)
        
        self.diff+=np.outer([neg(x) for x in vector],self.diff[end,:]) # subtract all precompositions with the differential (rows of diff)
        self.diff[:,end]+=np.dot(self.diff,vector) # add all postcompositions with the differential (columns of diff)
    
    def isolate_arrow(self,start, end, alg):
        """ Try to isotope away all arrows with the same start or the same end as 'arrow'.
        'arrow' is a list [start, end, alg], where 'alg' is a BNmor with a single entry.

        """
                
        face=alg.pairs[0][0]
        inverse_coeff=inverse(alg.pairs[0][1],self.field)
        
        def find_isotopy(index,entry):# Note: We are assuming here that there is at most one label to remove. 
            for pair in entry.pairs:
                if (pair[0]*face>0):
                    if ((self.diff)[index,end] is 0) & (index != end):
                        return BNmor([[pair[0]-face,pair[1]*inverse_coeff]],self.field)
            return 0 # zeor algebra element
        
        # first remove all other arrows with the same start
        for index in range(len(self.gens)):
            if ((self.diff)[index,end] is 0) and (index != end): # check that isotopy is valid.
                if self.diff[index,start] is not 0:
                    for pair in self.diff[index,start].pairs:
                        if pair[0]*face>0:# same face
                            self.isotopy(end,index,BNmor([[pair[0]-face,pair[1]*inverse_coeff]],self.field),"unsafe")
        # attempt to speed up , but actually MUCH slower...
        #isotopy_vector=np.array([find_isotopy(index,entry) for index,entry in enumerate(self.diff[:,start])])
        #self.isotopy_via_vector_end(end,isotopy_vector)
        
        # secondly, remove all remaining arrows with the same end
        for index in range(len(self.gens)):
            if ((self.diff)[start,index] is 0) & (index != start): # check that isotopy is valid.
                if self.diff[end,index] is not 0:
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
                    if self.diff[end,start_current] is not 0:
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
                    if self.diff[end_current,start] is not 0:
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
        matrix_D=np.array([[((entry is not 0) and entry.contains_D()) for entry in row] for row in self.diff])
        matrix_S=np.array([[((entry is not 0) and entry.contains_S()) for entry in row] for row in self.diff])
        size=len(matrix_D)
        return  all([list(matrix_D[:,index]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_D[index,:]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_S[:,index]).count(True)<2 for index in range(size)]) & \
                all([list(matrix_S[index,:]).count(True)<2 for index in range(size)])
    
    def to_multicurve(self):
        """convert a loop-type BNComplex into multicurve, ie split along differentials
        """
        if self.is_looptype == False:
            raise Exception('I cannot convert this complex into a multicurve, because it is not loop-type!')
        
        curves=[]
        index_remaining = list(range(len(self.gens)))
        
        while index_remaining != []:
            #split off a single curve
            curve = []
            current = index_remaining[0]
            
            while (current not in curve):
                index_remaining.remove(current)
                curve.append(current)
                try:# try to find an adjacent generator that is not already in 'curve'
                    # if unsuccessful, exit while loop 
                    def non_zero(index):
                        return (self.diff[index,current] is not 0) or (self.diff[current,index] is not 0)
                    index=find_first(index_remaining,non_zero)
                    current=index     
                except: 
                    pass
            current = curve[0] # go back to the start of the curve, in case the curve is an arc
            
            try:# try to find an adjacent generator that is not already in 'curve'
                # if unsuccessful, we have found the curve
                index=find_first(index_remaining,non_zero)
                current=index 
                while (current not in curve) or (current == curve[0]):
                    index_remaining.remove(current)
                    curve.insert(0,current)
                    try:# try to find an adjacent generator that is not already in 'curve'
                        # if unsuccessful, exit while loop, we have found the curve
                        index=find_first(index_remaining,non_zero)
                        current=index 
                    except: 
                        pass
            except: 
                pass
            curves.append(curve)
        
        genss=[[self.gens[index] for index in curve] for curve in curves]
        diffs=[[[self.diff[j,i] for j in curve] for i in curve] for curve in curves]
        return multicurve([BNComplex(gens,diff,self.field) for gens,diff in zip(genss,diffs)])
        
    
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
            print("iteration: "+str(iter)+" ("+str(round(time2-time1,1))+" sec)", end='\r')# testing how to monitor a process
            time1=time2
            #time.sleep(1)
        else:
            print("Clean-up: Terminated as a result of reaching maximum iteration value of", max_iter)
        
    def shift_h(self,shift): # create new complex
        new_gens = [gen.shift_h(shift) for gen in self.gens]
        if shift % 2 == 0:
            new_diff = self.diff
        elif shift % 2 == 1:
            def neg(alg):
                if alg is 0:
                    return 0
                return alg.negative(self.field)
            new_diff = [[neg(alg) for alg in row] for row in self.diff]
        else:
            raise Exception('Why are you trying to shift homological grading by something other than an integer? I cannot do that!')
        return BNComplex(new_gens,new_diff,self.field)
    
    def shift_q(self,shift): # modify the complex
        self.gens = [gen.shift_q(shift) for gen in self.gens]
    
    def cone(self,Hpower):
        
        shifted_complex = self.shift_h(1)
        shifted_complex.shift_q(2*Hpower)
        
        zero_matrix=np.array([[0 for i in range(len(self.gens))] for i in range(len(self.gens))])
        
        def fill_diagonal(i,j):
            if i==j:
                return BNmor([[Hpower,1],[-2*Hpower,(-1)**Hpower]],self.field)
            else:
                return 0
        
        Hdiagonal=np.array([[fill_diagonal(i,j) for j in range(len(self.gens))] for i in range(len(self.gens))])
        
        new_diff = np.concatenate(\
                (np.concatenate((self.diff,zero_matrix),axis=1),\
                np.concatenate((Hdiagonal,shifted_complex.diff),axis=1)),axis=0)
        
        return BNComplex(self.gens+shifted_complex.gens,new_diff,self.field)

def CobordismToBNAlg(cob,field=2):
    """ Convert a cobordism between two (1,3)-tangles into an element of BNAlgebra."""
    
    if cob==0:
        return 0# zero algebra element
    
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
        return BNobj(0,clt.q,clt.h) #b
    elif clt.arcs[0]==1:
        return BNobj(1,clt.q,clt.h) #c

def CobComplex2BNComplex(complex,field=2):
    gens=[CLT2BNObj(clt) for clt in complex.gens]
    diff=[[CobordismToBNAlg(cob,field) for cob in row] for row in complex.diff]
    BNcx = BNComplex(gens,diff,field)
    BNcx.eliminateAll()
    return BNcx

def BNObj2CLT(bnobj):
    arcs = []
    if bnobj.idem == 0:
        arcs = [3,2,1,0]
    else:
        arcs = [1,0,3,2]
    return CLT(1,3, arcs, bnobj.h, bnobj.q, bnobj.delta)

def BNAlg2Cob(morphism, sourceCLT, targetCLT):
    decos = []
    if morphism is 0:
        return ZeroCob
    
    for pair in morphism.pairs:
        if pair[0] > 0 : 
            decos.append([pair[0]-1, 0, 1, pair[1]])
        elif pair[0] == 0: 
            decos.append([0, 0, 0, pair[1]])
        elif pair[0] %2 == 0:
            power = int(Fraction(-1*pair[0], 2))
            decos.append([power, 0, 0, ((-1)** power) *pair[1]]) #(-H)^(n/2)
            decos.append([power - 1, 0, 1, ((-1)** (power -1)) *pair[1]]) #(-H)^(n/2-1) D
        elif pair[0] %2== 1: 
            power = int(Fraction(-1*pair[0] -1, 2))
            decos.append([power, 0, ((-1)** power)*pair[1] ])#(-H)^((n-1)/2) S
        else: 
            raise Exception("pair is not an integer?")
    Cob = Cobordism(sourceCLT, targetCLT, decos)
    # print("COB", Cob.front.top, Cob.front.bot)
    Cob.ReduceDecorations()
    return Cob
    
def BNComplex2CobComplex(BNcomplex):
    # BNcomplex.print()
    gens = [BNObj2CLT(bnobj) for bnobj in BNcomplex.gens]
    diff = [[BNAlg2Cob(morphism, gens[source], gens[target]) \
                 for source, morphism in enumerate(row)] for target, row in enumerate(BNcomplex.diff)]
    complex = CobComplex(gens, diff)
    # complex.print("old long")
    return complex
    
class multicurve(object):
    """A multicurve is a collection of loop-type BNcomplexes which are assumed to be connected.
    """
    def __init__(self,comps):
        self.comps=comps
    
    def draw(self, filename, vertex_switch="index_qhdelta",tangle=None):
        """create a pdf file which draws the graph of each component on a separate page; if tangle is specified, the first page is the tangle input.
        """
        for i,comp in enumerate(self.comps):
            comp.draw(filename+str(i)+".pdf",vertex_switch)
        #sleep(4)
        if tangle==None:
            tanglestr=""
        else:
            tanglestr="Output/"+filename+"_tangle.pdf "
            drawtangle(tangle,filename+"_tangle","slices",1)
        run("pdftk "+tanglestr+"".join(["Output/"+filename+str(i)+".pdf " for i in range(len(self.comps))])+"output Output/"+filename+".pdf", shell=True)
        run("rm "+tanglestr+"".join(["Output/"+filename+str(i)+".pdf " for i in range(len(self.comps))]), shell=True)

#todo: implement recognition of local systems (optional)
#todo: implement pairing theorem (just for fun!)

