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
import numpy as np
import pandas as pd
from fractions import Fraction
from random import choice
from subprocess import run
from tabulate import tabulate
from time import time

import Drawing
import BNAlgebra
import CobComplexes
from aux import *

class BNComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of BNAlgebra.obj as labels on the vertices.
        - A matrix of BNAlgebra.mor as the adjacency matrix.
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
            if mor == 0:
                return 0
            return mor.simplify_mor(field)
        self.diff = np.array([[simplify(mor) for mor in row] for row in self.diff])
    
    def __repr__(self):
        return "BNComplex({},{},{})".format(self.gens,[list(row) for row in self.diff],self.field)
    
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
        def ToStr(entry):
            if entry == 0:
                return ""
            return str(entry)
        string+=str(tabulate(pd.DataFrame([[ToStr(entry) for entry in row] for row in self.diff]),range(len(self.diff)),tablefmt="fancy_grid"))
        return string
    
    def save(self,filename):
        with open("examples/data/BNComplexes/"+filename, "w") as text_file:
            print(repr(self), file=text_file)

    def validate(self):
        # return True
        length = len(self.gens)
        if len(self.diff) != length:
            raise Exception('Differential does not have n rows (where n is the number of gens in chain complex)')
        for i in self.diff:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of gens in chain complex)')
        for i in range(length):
            if self.diff[i][i] != 0:
                raise Exception('Differential has self loops')
        
        squared = np.tensordot(self.diff,self.diff, axes=(-2,-1))
        for i,row in enumerate(squared):
            for j,mor in enumerate(row):
                if mor != 0:
                    mor.simplify_mor(self.field)
                    if mor != 0:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: Found non-zero term "+str(mor)+" in d² in row "+str(i)+" and column "+str(j)+".")
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Differential does not square to 0')
        
    def findIsom(self): 
        """Returns the location of the first isomorphism it finds
           If no isomorphism is found, returns None"""
        for targetindex, row in enumerate(self.diff):
            for sourceindex, morphism in enumerate(row):
                if (morphism != 0) and (morphism.is_isomorphism()):
                    return [sourceindex, targetindex]
        return None
    
    def eliminateIsom(self, sourceindex, targetindex):
        """Mutates self by eliminating the isomorphism at the specified location, via the gaussian elimination lemma
           Note that this does not check that the cobodism specified is actually an isomorphism.
        """
        Max=max(targetindex,sourceindex)
        Min=min(targetindex,sourceindex)
        label = self.diff[targetindex,sourceindex] # label of the cancelling arrow
        
        if (len(label.pairs)!=1) or (label.pairs[0][0]!=0):
            raise Exception('You cannot cancel this arrow!')
        
        coeff_inv=BNAlgebra.inverse(label.pairs[0][1],self.field) # inverse of the coefficient of the label of the cancelling arrow
        
        del self.gens[Max] # eliminate source and target from list of generators
        del self.gens[Min] # eliminate source and target from list of generators
        
        out_source = np.delete(self.diff[:,sourceindex],[Min,Max],0) #arrows starting at the source, omiting indices targetindex and sourceindex
        in_target = np.delete(self.diff[targetindex],[Min,Max],0) #arrows ending at the target, omiting indices targetindex and sourceindex
        in_target = np.array([(-entry)*coeff_inv for entry in in_target])
        
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
        self.validate()
       
    def isotopy(self,start,end,alg,switch="safe"):
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        
        if (switch=="safe") and ((self.diff)[start,end] != 0):
            raise Exception('This isotopy probably does not preserve the chain isomorphism type. There is an arrow going in the opposite direction of the isotopy.')
        self.diff[end,:]+=[(-alg)*element for element in self.diff[start,:]] # subtract all precompositions with the differential (rows of diff)
        self.diff[:,start]+=[element*alg for element in self.diff[:,end]] # add all postcompositions with the differential (columns of diff)
    
    def isotopy_via_vector_end(self,end,vector): # unused function: this is actually *much* slower
        """ Apply an isotopy along an arrow (start--->end) labelled by 'alg'.
        """
        self.diff+=np.outer([-x for x in vector],self.diff[end,:]) # subtract all precompositions with the differential (rows of diff)
        self.diff[:,end]+=np.dot(self.diff,vector) # add all postcompositions with the differential (columns of diff)
    
    def isolate_arrow(self,start, end, alg):
        """ Try to isotope away all arrows with the same start or the same end as 'arrow'.
        'arrow' is a list [start, end, alg], where 'alg' is a BNAlgebra.mor with a single entry.

        """
                
        face=alg.pairs[0][0]
        inverse_coeff=BNAlgebra.inverse(alg.pairs[0][1],self.field)
        
        def find_isotopy(index,entry):# Note: We are assuming here that there is at most one label to remove. 
            for pair in entry.pairs:
                if (pair[0]*face>0):
                    if ((self.diff)[index,end] == 0) & (index != end):
                        return BNAlgebra.mor([[pair[0]-face,pair[1]*inverse_coeff]],self.field)
            return 0 # zeor algebra element
        
        # first remove all other arrows with the same start
        for index in range(len(self.gens)):
            if ((self.diff)[index,end] == 0) and (index != end): # check that isotopy is valid.
                if self.diff[index,start] != 0:
                    for pair in self.diff[index,start].pairs:
                        if pair[0]*face>0:# same face
                            self.isotopy(end,index,BNAlgebra.mor([[pair[0]-face,pair[1]*inverse_coeff]],self.field),"unsafe")
        # attempt to speed up , but actually MUCH slower...
        #isotopy_vector=np.array([find_isotopy(index,entry) for index,entry in enumerate(self.diff[:,start])])
        #self.isotopy_via_vector_end(end,isotopy_vector)
        
        # secondly, remove all remaining arrows with the same end
        for index in range(len(self.gens)):
            if ((self.diff)[start,index] == 0) & (index != start): # check that isotopy is valid.
                if self.diff[end,index] != 0:
                    for pair in self.diff[end,index].pairs:
                        if pair[0]*face>0:# same face
                            self.isotopy(index,start,BNAlgebra.mor([[pair[0]-face,(-1)*pair[1]*inverse_coeff]],self.field),"unsafe")
    
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
                    if self.diff[end,start_current] != 0:
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
                    if self.diff[end_current,start] != 0:
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
                self.isolate_arrow(start_current, end_current, BNAlgebra.mor([self.diff[end_current,start_current].pairs[index_current]],self.field))
                #print("start:", start_current,"end:", end_current, "isotopy:",self.diff[end_current,start_current].pairs[index_current])
                remaining.remove(start_current)
                remaining.remove(end_current)
    
    def is_looptype(self):# todo: This could be made much more efficient!!
        """Return 'True' is the complex is loop-type.
        """
        matrix_D=np.array([[((entry != 0) and entry.contains_D()) for entry in row] for row in self.diff])
        matrix_S=np.array([[((entry != 0) and entry.contains_S()) for entry in row] for row in self.diff])
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
                        return (self.diff[index,current] != 0) or (self.diff[current,index] != 0)
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
        
    def clean_up(self,max_iter=1000):
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
            new_diff = [[-alg for alg in row] for row in self.diff]
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
                return BNAlgebra.mor([[Hpower,1],[-2*Hpower,(-1)**Hpower]],self.field)
            else:
                return 0
        
        Hdiagonal=np.array([[fill_diagonal(i,j) for j in range(len(self.gens))] for i in range(len(self.gens))])
        
        new_diff = np.concatenate(\
                (np.concatenate((self.diff,zero_matrix),axis=1),\
                np.concatenate((Hdiagonal,shifted_complex.diff),axis=1)),axis=0)
        
        return BNComplex(self.gens+shifted_complex.gens,new_diff,self.field)
    
    def ToCob(self):# does not work yet, since Z is not implemented over BNAlgebra and Cob is only implemented over Z.
        if self.field != 1:
            raise Exception("You are converting a complex over BNalgebra with coefficients in field={} into a complex over Cob. However, the category Cob is only implemented over integers, ie for field=1.".format(field))
        
        gens = [obj.ToCob() for obj in BNcomplex.gens]
        def convert(mor,source,target):
            if mor == 0:
                return 0
            return mor.ToCob(source,target)
        diff = [[convert(mor, gens[source], gens[target]) for source, mor in enumerate(row)] for target, row in enumerate(BNcomplex.diff)]
        complex = CobComplexes.CobComplex(gens, diff)
        # complex.print("old long")
        complex.validate()
        return complex
    
def importBNcx(filename):
    with open("examples/data/BNComplexes/"+filename, "r") as text_file:
        data = text_file.read()
        return eval(data)
    
class multicurve(object):
    """A multicurve is a collection of loop-type BNcomplexes which are assumed to be connected.
    """
    def __init__(self,comps):
        self.comps=comps
    
    def __repr__(self):
        return "multicurve({})".format(self.comps)
    
    def draw(self, filename, vertex_switch="qhdelta",title=["",""],tangle=None,thumbnails=False):
        """create a pdf file which draws the graph of each component on a separate page; if tangle is specified, the first page is the tangle input.
        """
        
        content="""\\documentclass{article}
\\usepackage[crop=off]{auto-pst-pdf}
\\usepackage{pst-node,rotating}
\\begin{document}
\\centering 
"""
        if thumbnails == True:
            vertex_switch=""
        
        for i,comp in enumerate(self.comps):
            size = len(comp.gens)
            if thumbnails == True:
                scaledsize=size/6
            else:
                scaledsize=size/3
            canvassize=scaledsize+0.5
            if "q" in vertex_switch:
                canvassize+=0.75
            if "h" in vertex_switch:
                canvassize+=0.75
            if "delta" in vertex_switch:
                canvassize+=0.75
            
            content+="\\begin{{pspicture}}(-{0},-{0})({0},{0})\n\\psset{{nodesep=2pt,shortput=nab,linewidth=1.5pt}}\n".format(canvassize)
        
            for i, gen in enumerate(comp.gens):
                if gen.idem == 0:
                    dottype="*" #black
                else:
                    dottype="" #white
                content+="\\cnode"+dottype+"("+str(scaledsize)+";"+str(360*i/size-90)+"){3.5pt}{P"+str(i)+"}\n"
                content+="\\nput{"+str(360*i/size-90)+"}{P"+str(i)+"}{$"+gen.grading2TeX(vertex_switch)+"$}\n"
        
            for i, diff_row in enumerate(comp.diff):
                for j, diff in enumerate(diff_row):
                    if diff !=0:
                        label=diff.label2TeX()
                        if thumbnails==True:
                            if ("D" in label) and ("S" in label):
                                content+="\\ncarc[arcangle=40,linestyle=dotted,dotsep=1pt]{P"+str(j)+"}{P"+str(i)+"}\n"
                                content+="\\ncarc[arcangle=-40]{P"+str(j)+"}{P"+str(i)+"}\n"
                            elif "D" in label:
                                content+="\\ncline[arrows=->]{P"+str(j)+"}{P"+str(i)+"}\n"
                            elif "S" in label:
                                content+="\\ncline[arrows=->,linestyle=dotted,dotsep=1pt]{P"+str(j)+"}{P"+str(i)+"}\n"
                        else:
                            if ("D" in label) and ("S" in label):
                                content+="\\ncarc[arcangle=40,arrows=->,linestyle=dotted,dotsep=1pt]{P"+str(j)+"}{P"+str(i)+"}\\naput{$"+diff.label2TeX("S")+"$}\n"
                                content+="\\ncarc[arcangle=-40,arrows=->]{P"+str(j)+"}{P"+str(i)+"}\\nbput{$"+diff.label2TeX("D")+"$}\n"
                            elif "D" in label:
                                content+="\\ncline[arrows=->]{P"+str(j)+"}{P"+str(i)+"}\\ncput*{$"+label+"$}\n"
                            elif "S" in label:
                                content+="\\ncline[arrows=->,linestyle=dotted,dotsep=1pt]{P"+str(j)+"}{P"+str(i)+"}\\ncput*{$"+label+"$}\n"
        
            content+="\\end{pspicture}\n"
        
        content+="\\end{document}"

        with open("examples/PSTricks/"+filename+".tex", "w") as text_file:
            print(content, file=text_file)
        
        run("cd examples/PSTricks && pdflatex -shell-escape '"+filename+".tex' > '"+filename+".out' 2>&1", shell=True)
        run("cd examples/PSTricks && rm "+(" ".join(["'"+filename+"."+string+"' " for string in ["log","aux","pdf","out"]])), shell=True)
        #subtitle="field="+str(comp.field)
        
        if tangle==None:
            tanglestr=""
        else:
            tanglestr="'examples/PSTricks/"+filename+"-tangle-pics.pdf'"
            if thumbnails:
                Drawing.drawtangle(tangle,filename,"plain",1)
            else:
                Drawing.drawtangle(tangle,filename,"slices",1,title)
        run("pdftk "+tanglestr+" 'examples/PSTricks/"+str(filename)+"-pics.pdf' output 'examples/"+filename+title[1]+".pdf'", shell=True)
        
#todo: implement recognition of local systems (optional)
#todo: implement pairing theorem (just for fun!)

