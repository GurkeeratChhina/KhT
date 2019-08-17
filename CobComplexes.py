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

import numpy as np
import pandas as pd
from tabulate import tabulate
from time import time

import Cobordisms as Cob
import Tangles

class CobComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of CLTS as labels on the vertices
        - A a matrix of cobordisms as the adjacency matrix.
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list 
        We assume that all entries of 'diff' are reduced in the sense that 'ReduceDecorations()' does not change them.
        """
    __slots__ = 'gens','diff'
    
    def __init__(self,gens,diff):
        self.gens = gens
        self.diff = np.array(diff)
    
    def print(self,switch="short"):
        """Print a complex in human readable form. The optional parameter should be one of the following strings: 
        - 'short' (default) prints only the length of cobordisms.
        - 'long' prints all cobordism data in a nice table.
        - 'old long' prints all cobordism data, but as a list of lists.
        """
        print("The generators:")
        print(pd.DataFrame({\
            "clt.pairs": [clt.pairs for clt in self.gens],\
            "q": [clt.q for clt in self.gens],\
            "h": [clt.h for clt in self.gens]
            },columns=["clt.pairs","q","h"]))
        print("The differential: ("+switch+" form)")
        #print(pd.DataFrame([[print(entry,switch)  for entry in row] for row in self.diff]))
        def prt(entry):
            if entry ==0:
                return ""
            return entry.print(switch)
        print(tabulate(pd.DataFrame([[prt(entry)  for entry in row] for row in self.diff]),range(len(self.diff)),tablefmt="fancy_grid"))
    
    def ValidMorphism(self): #checks that the differential squares to 0, has no self loops, and is a matrix of the correct size
        length = len(self.gens)
        if len(self.diff) != length:
            raise Exception('Differential does not have n rows (where n is the number of gens in chain complex)')
                
        for i in self.diff:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of gens in chain complex)')
        for i in range(length):
            if self.diff[i][i] is not 0:
                raise Exception('Differential has self loops')
        
        for i,row in enumerate(self.diff):
            for j,cob in enumerate(row):   
                if cob is not 0:
                    if cob.homogeneousQ() == False:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The component of the differential in row "+str(i)+" and column "+str(j)+" is not homoegenous:")
                        print(cob.print("long"))
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Non-homogeneous morphism in differential!')
                    if (self.gens[i]).h-(self.gens[j]).h!=1:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The homological grading along the component of the differential in row "+str(i)+" and column "+str(j)+" does not increase by 1:")
                        print(cob.print("long"))
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Something is wrong with the homological grading!')
                    if (self.gens[i]).q-(self.gens[j]).q+cob.deg()!=0:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The quantum grading is not preserved along the component of the differential in row "+str(i)+"(q:"+str((self.gens[i]).q)+") and column "+str(j)+"(q:"+str((self.gens[j]).q)+"):")
                        print(cob.print("long"), " degree:", cob.deg())
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Something is wrong with the quantum grading!')
        
        # Computing diff squared:
        squared = np.tensordot(self.diff,self.diff, axes=(-2,-1))
        for i,row in enumerate(squared):
            for j,cob in enumerate(row):
                if (cob is not 0) and (cob.ReduceDecorations() is not []):
                    print("!!!!!!!!!!!!!!!!!!")
                    print("ERROR: Found non-zero term in d² in row "+str(i)+" and column "+str(j)+":")
                    print(cob.print("long"))
                    print("!!!!!!!!!!!!!!!!!!")
                    raise Exception('Differential does not square to 0')
    
    def findIsom(self): 
        """Returns the location of the first isomorphism it finds
           If no isomorphism is found, returns None"""
        for targetindex, row in enumerate(self.diff):
            for sourceindex, cob in enumerate(row):
                if (cob != 0) and (cob.isIsom()):
                    return [sourceindex, targetindex]
        return None
    
    def eliminateIsom(self, sourceindex, targetindex):
        """Mutates self by eliminating the isomorphism at the specified location, via the gaussian elimination lemma
           Note that this does not check that the cobodism specified is actually an isomorphism.
        """
        
        Max=max(targetindex,sourceindex)
        Min=min(targetindex,sourceindex)
        
        del self.gens[Max] # eliminate source and target from list of generators
        del self.gens[Min] # eliminate source and target from list of generators
        
        out_source = np.delete(self.diff[:,sourceindex],[Min,Max],0) #arrows starting at the source, omiting indices targetindex and sourceindex
        in_target = np.delete(self.diff[targetindex],[Min,Max],0) #arrows ending at the target, omiting indices targetindex and sourceindex
        
        if (self.diff[targetindex,sourceindex]).decos[0][-1]==1: # add minus sign; in the case the coefficient is -1, the signs cancel.
            in_target=np.array([-entry for entry in in_target])
        
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
    
    def shift_qhd(self,q,h,delta):
        self.gens=[clt.shift_qhd(q,h,delta) for clt in self.gens]

def AddCapToCLT(clt, i, grshift = "false"): 
    """creates a new CLT which is clt with a cap added to it at index i
       Here 0 <= i <= tangle.bot
       If grshift is set to "true", then it also shifts the homological and quantum grading of the clt by 1 each
       This should only be used when adding a crossing"""
    def incrementby2(j): #increment TEI by 2 if greater than where cap is to be inserted
        if j >= clt.top+i:
            return j+2
        else:
            return j
    newarcs = [incrementby2(x) for x in clt.arcs[:clt.top+i]] + [clt.top+i +1, clt.top+i] + [incrementby2(x) for x in clt.arcs[clt.top+i:]]
    if grshift == "true": # shift grading
        newh = clt.h +1
        newq = clt.q +1
        newdelta = clt.delta -0.5
    else: # dont shift grading
        newh = clt.h
        newq = clt.q
        newdelta = clt.delta
    return Tangles.CLT(clt.top, clt.bot+2, newarcs, newh, newq, newdelta)

def AddCap(Complex, i, grshift = "false"):
    """ Creates a new complex by adding a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot
        If grshift is set to "true", then it applies a grading shift to every tangle (as above)
        Furthermore, it will flip the sign on every cobordism, which is the convention used so that the differential will square to 0 when adding a crossing"""
    Newgens = [AddCapToCLT(clt, i, grshift) for clt in Complex.gens]
    length = len(Complex.gens)
    Newdiff = []
    for target ,row in enumerate(Complex.diff):
        NewRow=[]
        for source, cob in enumerate(row):
            if cob == 0: #adding a cap to the zero cobordism does nothing
                NewRow.append(0)
            else: #not the zero cobordism
                def incrementby2(j):
                    if  j >= cob.front.top+i: #increment TEI by 2 if greater than where cap is to be inserted
                        return j+2
                    else:
                        return j
                newcomps=[[incrementby2(x) for x in comp] for comp in cob.comps] + [[cob.front.top +i, cob.front.top+i+1]] #shifts TEI of old components, and adds new component to the end of the list
                if grshift == "true":
                    newDecos = [ NewDeco[:-1] + [0] + [NewDeco[-1]*-1] for NewDeco in cob.decos] #adds the new component without a dot and flips the sign on the coefficient
                else: 
                    newDecos = [ NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos] #adds the new component without a dot
                NewCob = Cob.Cobordism(Newgens[source], Newgens[target], Cob.simplify_decos(newDecos), newcomps)
                NewRow.append(NewCob.ReduceDecorations())
        Newdiff.append(NewRow)
    return CobComplex(Newgens, Newdiff)

def AddCupToCLT(clt, i):
    """ Adds a cup to the clt at index i, where 0 <= i <= clt.bot -2
        Returns a list of 2 clt if there is a closed component, obtained by neckcutting
        Otherwise returns a list of 1 clt """
    newgens = []
    def decrementby2(j):#decrement TEI by 2 if greater than where cup is to be inserted
            if j >= clt.top +i:
                return j-2
            else:
                return j
    if clt.arcs[clt.top +i] == clt.top+i+1: # adding the cup makes a closed component
        newarcs = [decrementby2(x) for x in clt.arcs if x != clt.top+i and x!= clt.top+i+1] #removes the closed component from the arcs, shifts remaining TEI
        newgens.append(Tangles.CLT(clt.top, clt.bot -2, newarcs, clt.h, clt.q+1, clt.delta+0.5)) #neckcutting
        newgens.append(Tangles.CLT(clt.top, clt.bot -2, newarcs, clt.h, clt.q-1, clt.delta-0.5)) #neckcutting
    else: # adding the cup doesnt make closed components
        leftend = clt.arcs[clt.top + i] #the endpoint of the arc which connects at i
        rightend = clt.arcs[clt.top+i+1] #the endpoint of the arc which connects at i+1
        newarcs = clt.arcs.copy()
        newarcs[leftend] = rightend # connecting the two arcs via the cup
        newarcs[rightend] = leftend
        newarcs1 = [decrementby2(entry) for x, entry in enumerate(newarcs) if x != clt.top+i and x != clt.top+i+1] #removes the ends which don't exist anymore and shifts remaining TEI
        newgens.append(Tangles.CLT(clt.top, clt.bot -2, newarcs1, clt.h, clt.q, clt.delta))
    return newgens

def AddCup(Complex, i): # TODO: reduce decorations
    """ Here 0 <= i <= tangle.bot -2"""
    newgens = []
    for clt in Complex.gens:
        newgens.extend(AddCupToCLT(clt, i))
    Newdiff = []
    for target, row in enumerate(Complex.diff):
        newRow = []
        nextRow = []
        for source, cob in enumerate(row):
            def decrementby2(j): #shifts TEI by -2 if it is greater than where the cup i added
                if j >= cob.front.top+i:
                    return j-2
                else:
                    return j
            
            def compcup(component): #computes the new component after adding a cup, for components except those with only 2 entries, which are clt.top+i and clt.top+i+1
                return [decrementby2(j) for j in component if j != cob.front.top+i and j!= cob.front.top+i+1]
            def incrementH(decoration): # Increments H of a decoration by 1
                return [decoration[0]+1] + decoration[1:]
            def negativeH(decoration): # Increments H of a decoration by 1, and flips the sign on the coefficient
                return [decoration[0]+1] + decoration[1:-1] +[decoration[-1]*-1]
            def adddotonindex(decoration, index_to_add): # adds dot to the component labeled by index_to_add if there is no dot already, otherwise does nothing
                return decoration[:index_to_add+1] + [1] + decoration[index_to_add+2:]
            
            if Complex.gens[source].arcs[Complex.gens[source].top +i] == Complex.gens[source].top+i+1 \
                and Complex.gens[target].arcs[Complex.gens[target].top +i] == Complex.gens[target].top+i+1: # source and target are both closed, add 2 cobordisms to newRow and nextRow
                if cob == 0:
                    newRow.append(0)
                    newRow.append(0)
                    nextRow.append(0)
                    nextRow.append(0)
                else: # is not the 0 cob
                    magic_index = 0 # the index of the first component that the cup connects to
                    for x,comp in enumerate(cob.comps):
                        if cob.front.top +i in comp:
                            if len(comp) != 2:
                                raise Exception("closed component doesn't have 2 TEI")
                            magic_index = x
                            break
                    else:
                        raise Exception("no magic_index")
                    # magic_index = next(j for j, v in enumerate(cob.comps) if len(v) == 2 and cob.front.top+i in v and cob.front.top+i+1 in v) #computes index of closed component, should be the same as magic_index computed generally above
                    newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps if cob.front.top+i not in comp ]
                    def delclosedcomp(decoration): # deletes the decoration corresponding to the closed component, which is to be removed
                        return decoration[:magic_index+1] + decoration[magic_index+2:]
                    def computedeco4(decoration):
                        if decoration[magic_index +1] == 0:
                            return delclosedcomp(decoration) # id on all other components if no dot
                        else:
                            return delclosedcomp(incrementH(decoration)) # H times rest of cobordism if dot
                    
                    newDecos1 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index+1] == 0] # identity on all other components if no dot, otherwise 0
                    newDecos3 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index+1] == 1] # identity on all other components if dot, otherwise 0
                    newDecos4 = [computedeco4(deco) for deco in cob.decos] # see computedeco4()
                    newsourceclts = AddCupToCLT(cob.front, i)
                    newtargetclts = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cob.Cobordism(newsourceclts[0], newtargetclts[0], Cob.simplify_decos(newDecos1), newcomps)
                    Cobordism2 = 0 # The neckcutting/delooping isomorphism will give us 0 for Cobordism 2
                    Cobordism3 = Cob.Cobordism(newsourceclts[0], newtargetclts[1], Cob.simplify_decos(newDecos3), newcomps)
                    Cobordism4 = Cob.Cobordism(newsourceclts[1], newtargetclts[1], Cob.simplify_decos(newDecos4), newcomps)
                    newRow.append(Cobordism1.ReduceDecorations())
                    newRow.append(Cobordism2)
                    nextRow.append(Cobordism3.ReduceDecorations())
                    nextRow.append(Cobordism4.ReduceDecorations())
            elif Complex.gens[target].arcs[Complex.gens[target].top +i] == Complex.gens[target].top+i+1: # source is open but target is closed, add 1 cobordism to each
                if cob == 0:
                    newRow.append(0)
                    nextRow.append(0)
                else: # is not the 0 cob
                    magic_index = 0 # the index of the first component that the cup connects to
                    for x,comp in enumerate(cob.comps):
                        if cob.front.top +i in comp:
                            magic_index = x
                            break
                    else:
                        raise Exception("no magic_index")
                    newcomps = [compcup(component) for component in cob.comps]
                    newDecos1 = []
                    for decoration in cob.decos:
                        if decoration[magic_index+1] == 0: # is D - H if no dot, 0 otherwise
                            newDecos1.extend([adddotonindex(decoration, magic_index), negativeH(decoration)])
                    newDecos2 = [deco for deco in cob.decos] # is identity regardless of dot or not
                    newsourceclt = AddCupToCLT(cob.front, i)[0]
                    newtargetclts = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cob.Cobordism(newsourceclt, newtargetclts[0], Cob.simplify_decos(newDecos1), newcomps)
                    Cobordism2 = Cob.Cobordism(newsourceclt, newtargetclts[1], Cob.simplify_decos(newDecos2), newcomps)
                    newRow.append(Cobordism1.ReduceDecorations())
                    nextRow.append(Cobordism2.ReduceDecorations())
            elif Complex.gens[source].arcs[Complex.gens[source].top +i] == Complex.gens[source].top+i+1: # source is closed but target is open, add 2 cobordisms to newRow only
                if cob == 0:
                    newRow.append(0)
                    newRow.append(0)
                else: # is not the 0 cob
                    magic_index = 0 # the index of the first component that the cup connects to
                    for x,comp in enumerate(cob.comps):
                        if cob.front.top +i in comp:
                            magic_index = x
                            break
                    else:
                        raise Exception("no magic_index")
                    newcomps = [compcup(component) for component in cob.comps]
                    def computedeco2(decoration):
                        if decoration[magic_index+1] == 1:
                            return incrementH(decoration) # H trading if already a dot
                        else:
                            return adddotonindex(decoration, magic_index) # otherwise add a dot
                    newDecos1 = [deco for deco in cob.decos] # always the identity
                    newDecos2 = [computedeco2(deco) for deco in cob.decos]
                    newsourceclts = AddCupToCLT(cob.front, i)
                    newtargetclt = AddCupToCLT(cob.back, i)[0]
                    Cobordism1 = Cob.Cobordism(newsourceclts[0], newtargetclt, Cob.simplify_decos(newDecos1), newcomps)
                    Cobordism2 = Cob.Cobordism(newsourceclts[1], newtargetclt, Cob.simplify_decos(newDecos2), newcomps)
                    newRow.append(Cobordism1.ReduceDecorations())
                    newRow.append(Cobordism2.ReduceDecorations())
            else: # source and target are both open, add 1 cobordism to newRow only
                if cob == 0:
                    newRow.append(0)
                else: # is not the 0 cob
                    # print("i", i)
                    # print("cob.front.top", cob.front.top)
                    # print("cob.front.top+i", cob.front.top+i)
                    # print(cob.comps)
                    magic_index = 0 # the index of the first component that the cup connects to
                    for x1,comp in enumerate(cob.comps):
                        if cob.front.top +i in comp:
                            magic_index = x1
                            break
                    else:
                        raise Exception("no magic_index")
                    magic_index_2 = 0 # index of component containing i+1
                    for x2,comp in enumerate(cob.comps):
                        if cob.front.top+i+1 in comp:
                            magic_index_2 = x2
                            break
                    else:
                        raise Exception("no magic_index_2")
                    # print("magic indices:", magic_index, magic_index_2)
                    newDecos1 = []
                    if magic_index == magic_index_2: # only one component being connected via cup
                        # print("magic indices are the same")
                        comp = cob.comps[magic_index] # the component being connected
                        x1 = comp.index(cob.front.top +i) 
                        x2 = comp.index(cob.front.top +i + 1)
                        comp1 = comp[:min(x1, x2)] + comp[max(x1, x2)+1:]
                        comp2 = comp[min(x1, x2) +1:max(x1, x2)]
                        newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps[:magic_index] +[comp1, comp2] + cob.comps[magic_index+1:]]
                        def computedeco1comps(decoration):
                            if decoration[magic_index+1] == 1:
                                return [decoration[:magic_index+1] +[1] + decoration[magic_index+1:]]
                            else:
                                return [decoration[:magic_index+1] +[1] + decoration[magic_index+1:] , \
                                        decoration[:magic_index+2] +[1] + decoration[magic_index+2:], \
                                        [decoration[0]+1] + decoration[1:magic_index+1] +[0] + decoration[magic_index+1:-1] + [decoration[-1]*-1]]
                        for deco in cob.decos:
                            newDecos1.extend(computedeco1comps(deco))
                    else: # two seperate components being connected via cup
                        comp1 = cob.comps[magic_index] # the component containing i
                        comp2 = cob.comps[magic_index_2] # the component containing i+1
                        location1 = comp1.index(cob.front.top+i)
                        location2 = comp2.index(cob.front.top+i+1)
                        comp2 = (comp2[location2+1:] + comp2[:location2]) #rotates list until i+1 is at front, and removes it
                        if location1 %2 == location2%2: # if top/bot dont line up, flip comp2
                            comp2.reverse()
                        comp1 = comp1[:location1] + comp2 + comp1[location1+1:] # insert comp2 into comp1, and dont include element i
                        newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps[:magic_index] + [comp1] + cob.comps[magic_index+1:]]
                        del newcomps[magic_index_2]
                        def computedeco2comps(decoration):
                            if decoration[magic_index+1] == 1 and decoration[magic_index_2+1]==1:
                                # newdecos = incrementH(decoration)
                                # del newdecos[magic_index_2+1]
                                # print("old decos:", incrementH(decoration)[:magic_index_2+1] + incrementH(decoration)[magic_index_2+2:])
                                # print("newdecos:", newdecos)
                                return incrementH(decoration)[:magic_index_2+1] + incrementH(decoration)[magic_index_2+2:]
                            elif decoration[magic_index+1] == 1 or decoration[magic_index_2+1]==1:
                                # newdecos = adddotonindex(decoration, magic_index)
                                # del newdecos[magic_index_2+1]
                                # print("old decos:", adddotonindex(decoration, magic_index)[:magic_index_2+1] + adddotonindex(decoration, magic_index)[magic_index_2+2:])
                                # print("newdecos:", newdecos)
                                return adddotonindex(decoration, magic_index)[:magic_index_2+1] + adddotonindex(decoration, magic_index)[magic_index_2+2:]
                            else:
                                return decoration[:magic_index_2+1] + decoration[magic_index_2+2:]
                        newDecos1 = [computedeco2comps(deco) for deco in cob.decos]
                    newsourceclt = AddCupToCLT(cob.front, i)[0]
                    newtargetclt = AddCupToCLT(cob.back, i)[0]
                    Cobordism1 = Cob.Cobordism(newsourceclt, newtargetclt, Cob.simplify_decos(newDecos1), newcomps)
                    newRow.append(Cobordism1.ReduceDecorations())
        Newdiff.append(newRow)
        if Complex.gens[target].arcs[Complex.gens[target].top +i] == Complex.gens[target].top+i+1: # target is closed
            Newdiff.append(nextRow)
    return CobComplex(newgens, Newdiff)
  
def AddPosCrossing(Complex, i):
    CapCup = AddCap(AddCup(Complex, i), i, "true")
    sourcegens = Complex.gens
    targetgens = CapCup.gens
    Newgens = sourcegens + targetgens
    TopLeft = Complex.diff
    BottomRight = CapCup.diff
    length1 = len(Complex.gens)
    length2 = len(CapCup.gens)
    TopRight = np.full((length1, length2), 0, Cob.Cobordism)
    BottomLeft = []
    for x, targetclt in enumerate(Complex.gens):
        if targetclt.arcs[targetclt.top +i] == targetclt.top+i+1: #targetclt is closed
            newRow = []
            nextRow = []
            for y, sourceclt in enumerate(Complex.gens):
                if x == y:
                    newTarget1 = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
                    newTarget2 = AddCapToCLT(AddCupToCLT(targetclt, i)[1], i, "true")
                    newcomps = Tangles.components(sourceclt, newTarget1)
                    magic_index = 0
                    for z,comp in enumerate(newcomps):
                        if sourceclt.top +i in comp:
                            magic_index = z
                            break
                    decos1 = [[0] + [0 for comp in newcomps[:magic_index]] + [1] + [0 for comp in newcomps[magic_index+1:]] + [1], [1] + [0 for comp in newcomps] + [-1]] #Compute new decos via neckcutting
                    decos2 = [[0] + [0 for comp in newcomps] + [1]]
                    NewCobordism1 = Cob.Cobordism(sourceclt, newTarget1, decos1, newcomps)
                    NewCobordism2 = Cob.Cobordism(sourceclt, newTarget2, decos2, newcomps)
                    newRow.append(NewCobordism1)
                    nextRow.append(NewCobordism2)
                else:
                    newRow.append(0)
                    nextRow.append(0)
            BottomLeft.append(newRow)
            BottomLeft.append(nextRow)
        else:
            newRow = []
            for y, sourceclt in enumerate(Complex.gens):
                if x == y:
                    newTarget = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
                    decos = [[0] + [0 for comp in Tangles.components(sourceclt, newTarget)] + [1]]
                    NewCobordism = Cob.Cobordism(sourceclt, newTarget, decos)
                    newRow.append(NewCobordism)
                else:
                    newRow.append(0)
            BottomLeft.append(newRow)
    Newdiff = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewComplex = CobComplex(Newgens, Newdiff)
    return NewComplex

def grshiftclt(clt):
    return Tangles.CLT(clt.top, clt.bot, clt.arcs, clt.h+1, clt.q+1, clt.delta-0.5)

def grshiftcob(cob):
    if cob == 0:
        return 0
    newDecos = [deco[:-1] + [deco[-1]*-1] for deco in cob.decos]
    return Cob.Cobordism(grshiftclt(cob.front), grshiftclt(cob.back), newDecos, cob.comps)
    
def AddNegCrossing(Complex, i):
    targetgens = [grshiftclt(clt) for clt in Complex.gens]
    CapCup = AddCap(AddCup(Complex, i), i)
    sourcegens = CapCup.gens
    Newgens = sourcegens + targetgens
    TopLeft = CapCup.diff
    BottomRight = [[grshiftcob(cob) for cob in row] for row in Complex.diff]
    length1 = len(sourcegens)
    length2 = len(targetgens)
    TopRight = np.full((length1, length2), 0, Cob.Cobordism)
    BottomLeft = [] 
    for x, targetclt in enumerate(Complex.gens):
        newRow = []
        for y, sourceclt in enumerate(Complex.gens):
            if sourceclt.arcs[sourceclt.top +i] == sourceclt.top+i+1: #sourceclt is closed
                if x == y:
                    newSource1 = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
                    newSource2 = AddCapToCLT(AddCupToCLT(sourceclt, i)[1], i)
                    newTarget = targetgens[x]
                    newcomps = Tangles.components(newSource2, newTarget)
                    magic_index = 0
                    for z,comp in enumerate(newcomps):
                        if sourceclt.top +i in comp:
                            magic_index = z
                            break
                    decos1 = [[0] + [0 for comp in newcomps] + [1]]
                    decos2 = [[0] + [0 for comp in newcomps[:magic_index]] + [1] + [0 for comp in newcomps[magic_index+1:]] + [1]]
                    newCobordism1 = Cob.Cobordism(newSource1, newTarget, decos1, newcomps)
                    newCobordism2 = Cob.Cobordism(newSource2, newTarget, decos2, newcomps)
                    newRow.append(newCobordism1)
                    newRow.append(newCobordism2)
                else:
                    newRow.append(0)
                    newRow.append(0)
            else: #sourceclt is open
                if x == y:
                    newSource = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
                    newTarget = targetgens[x]
                    decos = [[0] + [0 for comp in Tangles.components(newSource, newTarget)] + [1]]
                    newCobordism = Cob.Cobordism(newSource, newTarget, decos)
                    newRow.append(newCobordism)
                else:
                    newRow.append(0)
        BottomLeft.append(newRow)
    Newdiff = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewComplex = CobComplex(Newgens, Newdiff)
    return NewComplex

def BNbracket(string,pos=0,neg=0,start=1,options="unsafe"):
    """compute the Bar-Natan bracket for tangle specified by 'string', which is a concatenation of words <type>+<index>, separated by '.' read from right to left, for each elementary tangle slice, read from top to bottom, where:
    <type> is equal to:
        'pos': positive crossing
        'neg': negative crossing
        'cup': cup
        'cap': cap
    <index> is the index at which the crossing, cap or cup sits. 
    'pos' and 'neg' are the numbers of positive and negative crossings. 
    The first optional parameter 'start' is an integer which specifies the number of tangle ends at the top.
    The second optional parameter 'options' is either 'unsafe' (default) or 'safe'. The latter performs some sanity checks, but slows down the computation.
    E.g. 'BNbracket('cup0pos0',2)' is a (2,0)-tangle which is decomposed as a positive crossing followed by a cap.
    """
    stringlist=[[word[0:3],int(word[3:])] for word in string.split('.')]
    stringlist.reverse()
    cx=CobComplex([Tangles.CLT(start,start,[start+i for i in range(start)]+[i for i in range(start)], 0,0,0)], [[0]])
    print("Computing the Bar-Natan bracket for the tangle\n\n"+string+"\n\n"+"with "+str(start)+" ends at the top, "+str(pos)+\
          " positive crossings, "+str(neg)+" negative crossings and "+str(len(stringlist))+" slices in total.")
          
    time0=time()
    time1=time0
    
    for i,word in enumerate(stringlist):
        
        time2=time()
        print("slice "+str(i)+"/"+str(len(stringlist))+": adding "+word[0]+" at index "+str(word[1])+" to tangle. ("+str(len(cx.gens))+" objects, "+str(round(time2-time1,1))+" sec)", end='\r')# monitor \n ->\r
        time1=time2
        
        if word[0]=="pos":
            cx=AddPosCrossing(cx, word[1])
            #print("before eliminateAll")
            #cx.print, "old long")
            if options=="safe": cx.ValidMorphism()
            cx.eliminateAll()
            # print("after eliminateAll")
            # cx.print, "old long")
        
        if word[0]=="neg":
            cx=AddNegCrossing(cx, word[1])
            #print("before eliminateAll")
            #cx.print, "old long")
            if options=="safe": cx.ValidMorphism()
            cx.eliminateAll()
            # print("after eliminateAll")
            # cx.print, "old long")
        
        if word[0]=="cup":
            cx=AddCup(cx, word[1])
            #print("before eliminateAll")
            #cx.print, "old long")
            if options=="safe": cx.ValidMorphism()
            cx.eliminateAll()
            #print("after eliminateAll")
            #cx.print, "old long")
        
        if word[0]=="cap":
            cx=AddCap(cx, word[1])
            if options=="safe": cx.ValidMorphism()
        
        

    cx.shift_qhd(pos-2*neg,-neg,0.5*neg)
    
    print("Completed the computation successfully after "+str(round(time1-time0,1))+" second(s).                        ")
    return cx


