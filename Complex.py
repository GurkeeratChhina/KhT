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
import math
from KhT import *
from Tangles import *
from Cobordisms import *
from Drawing import *
import time

class ChainComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of CLTS as labels on the vertices
        - A a matrix of cobordisms as the adjacency matrix.
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list """
    def __init__(self,listofclt,morphisms):
        self.elements = listofclt
        self.morphisms = morphisms
    
    def ValidMorphism(self): #checks that the differential squares to 0, has no self loops, and is a matrix of the correct size
        length = len(self.elements)
        if len(self.morphisms) != length:
            raise Exception('Differential does not have n rows (where n is the number of elements in chain complex)')
        for i in self.morphisms:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of elements in chain complex)')
        for i in range(length):
            if self.morphisms[i][i].decos != []:
                raise Exception('Differential has self loops')
        
        squared = np.tensordot(self.morphisms,self.morphisms, axes=(-2,-1))
        for i in flatten(squared):
            if i.ReduceDecorations() != []:
                raise Exception('Differential does not square to 0')
    
    def findIsom(self): 
        """Returns the location of the first isomorphism it finds
           If no isomorphism is found, returns None"""
        for targetindex, row in enumerate(self.morphisms):
            for sourceindex, cob in enumerate(row):
                if cob.isIsom():
                    return [sourceindex, targetindex]
        return None
    
    def eliminateIsom(self, sourceindex, targetindex):
        """Mutates self by eliminating the isomorphism at the specified location, via the gaussian elimination lemma
           Note that this does not check that the cobodism specified is actually an isomorphism"""
        Isom = self.morphisms[targetindex][sourceindex]
        Inverse = Cobordism(Isom.back, Isom.front, Isom.decos, Isom.comps) #because the isomorphism is +- id, taking the inverse (ie flipping the cobordism upside down) amounts to swapping the front and back
        newElements = [elem for x, elem in enumerate(self.elements) if x != sourceindex and x != targetindex] # both the source and the target are eliminated from the matrix via the lemma
        newMorphisms = []
        for target, row in enumerate(self.morphisms):
            newRow = []
            if target != sourceindex and target != targetindex: 
                for source, cob in enumerate (row):
                    if source != sourceindex and source != targetindex:
                        temp1 = (self.morphisms[targetindex][source]*Inverse)*self.morphisms[target][sourceindex] #composes the 3 morphisms as in the lemma
                        newCob = cob + temp1.negative()
                        newRow.append(newCob)
                newMorphisms.append(newRow)
        self.elements = newElements
        self.morphisms = newMorphisms

    def eliminateAll(self): #mutates self by eliminating isomorphisms as long as it can find one
        while True:
            index_to_eliminate = self.findIsom()
            if index_to_eliminate is None:
                break
            else:
                self.eliminateIsom(index_to_eliminate[0], index_to_eliminate[1])
    
    def shift_qhd(self,q,h,delta):
        self.elements=[clt.shift_qhd(q,h,delta) for clt in self.elements]

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
    if grshift == "true":
        newpgr = clt.pgr +1
        newqgr = clt.qgr +1
        newdgr = clt.dgr -0.5
    else:
        newpgr = clt.pgr
        newqgr = clt.qgr
        newdgr = clt.dgr
    return CLT(clt.top, clt.bot+2, newarcs, newpgr, newqgr, newdgr)

def AddCap(Complex, i, grshift = "false"):
    """ Creates a new complex by adding a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot
        If grshift is set to "true", then it applies a grading shift to every tangle (as above)
        Furthermore, it will flip the sign on every cobordism, which is the convention used so that the differential will square to 0 when adding a crossing"""
    NewElements = [AddCapToCLT(clt, i, grshift) for clt in Complex.elements]
    length = len(Complex.elements)
    NewMorphisms = []
    for j ,row in enumerate(Complex.morphisms):
        NewRow=[]
        for k, cob in enumerate(row):
            if cob.decos == []: #adding a cap to the zero cobordism does nothing
                NewRow.append(ZeroCob)
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
                NewCob = Cobordism(NewElements[k], NewElements[j], newDecos, newcomps)
                NewRow.append(NewCob)
        NewMorphisms.append(NewRow)
    return ChainComplex(NewElements, NewMorphisms)

def AddCupToCLT(clt, i):
    """ Adds a cup to the clt at index i, where 0 <= i <= clt.bot
        Returns a list of 2 clt if there is a closed component, obtained by neckcutting
        Otherwise returns a list of 1 clt """
    newElements = []
    def decrementby2(j):#decrement TEI by 2 if greater than where cup is to be inserted
            if j >= clt.top +i:
                return j-2
            else:
                return j
    if clt.arcs[clt.top +i] == clt.top+i+1: # adding the cup makes a closed component
        newarcs = [decrementby2(x) for x in clt.arcs if x != clt.top+i and x!= clt.top+i+1] #removes the closed component from the arcs, shifts remaining TEI
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, clt.pgr, clt.qgr+1, clt.dgr+0.5)) #neckcutting
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, clt.pgr, clt.qgr-1, clt.dgr-0.5)) #neckcutting
    else: # adding the cup doesnt make closed components
        leftend = clt.arcs[clt.top + i] #the endpoint of the arc which connects at i
        rightend = clt.arcs[clt.top+i+1] #the endpoint of the arc which connects at i+1
        newarcs = clt.arcs.copy()
        newarcs[leftend] = rightend # connecting the two arcs via the cup
        newarcs[rightend] = leftend
        newarcs1 = [decrementby2(entry) for x, entry in enumerate(newarcs) if x != clt.top+i and x != clt.top+i+1] #removes the ends which don't exist anymore and shifts remaining TEI
        newElements.append(CLT(clt.top, clt.bot -2, newarcs1, clt.pgr, clt.qgr, clt.dgr))
    return newElements

def AddCup(Complex, i): # TODO: reduce decorations
    """ Here 0 <= i <= tangle.bot -2"""
    newElements = []
    for clt in Complex.elements:
        newElements.extend(AddCupToCLT(clt, i))
    NewMorphisms = []
    for j, row in enumerate(Complex.morphisms): # j is the target clt
        newRow = []
        nextRow = []
        for k, cob in enumerate(row): # k is the source clt
            def decrementby2(j): #shifts TEI by -2 if it is greater than where the cup i added
                if j > cob.front.top+i+1:
                    return j-2
                else:
                    return j
            magic_index = 0
            if cob.decos != []:
                for x,comp in enumerate(cob.comps): #TODO: find magic index without iterating twice, dont need defining function
                    if cob.front.top +i in comp:
                        magic_index = x
                        break
            def compcup(component): #computes the new component after adding a cup, for components containing at least 3 ends
                return [decrementby2(j) for j in component if j != cob.front.top+i and j!= cob.front.top+i+1]
            def incrementH(decoration):
                return [decoration[0]+1] + decoration[1:]
            def negativeH(decoration):
                return [decoration[0]+1] + decoration[1:-1] +[decoration[-1]*-1]
            def adddotonmagicindex(decoration):
                return decoration[:magic_index+1] + [1] + decoration[magic_index+2:]
            
            if Complex.elements[k].arcs[Complex.elements[k].top +i] == Complex.elements[k].top+i+1 \
                and Complex.elements[j].arcs[Complex.elements[j].top +i] == Complex.elements[j].top+i+1: # source and target are both closed, add 2 cobordisms to newRow and nextRow
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    magic_index = cob.comps.index([cob.front.top+i, cob.front.top+i+1]) #computes index of closed component
                    newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps if cob.front.top+i not in comp ]
                    def delclosedcomp(decoration):
                        return decoration[:magic_index+1] + decoration[magic_index+2:]
                    def computedeco4(decoration):
                        if decoration[magic_index +1] == 0:
                            return delclosedcomp(decoration)
                        else:
                            return delclosedcomp(incrementH(decoration))
                    
                    newDecos1 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index+1] == 0]
                    newDecos3 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index+1] == 1]
                    newDecos4 = [computedeco4(deco) for deco in cob.decos]
                    NewDecos1 = [deco for deco in newDecos1 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    NewDecos3 = [deco for deco in newDecos3 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    NewDecos4 = [deco for deco in newDecos4 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    simplify_decos(NewDecos4)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], NewDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[1], bots[0], [], newcomps)
                    Cobordism3 = Cobordism(tops[0], bots[1], NewDecos3, newcomps)
                    Cobordism4 = Cobordism(tops[1], bots[1], NewDecos4, newcomps)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
                    nextRow.append(Cobordism3)
                    nextRow.append(Cobordism4)
            elif Complex.elements[j].arcs[Complex.elements[j].top +i] == Complex.elements[j].top+i+1: # source is open but target is closed, add 1 cobordism to each
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    newcomps = [compcup(component) for component in cob.comps]
                    newDecos1 = []
                    for decoration in cob.decos:
                        newDecos1.extend([adddotonmagicindex(decoration), negativeH(decoration)])
                    newDecos2 = [deco for deco in cob.decos]
                    NewDecos1 = [deco for deco in newDecos1 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    NewDecos2 = [deco for deco in newDecos2 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    simplify_decos(NewDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], NewDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[0], bots[1], NewDecos2, newcomps)
                    newRow.append(Cobordism1)
                    nextRow.append(Cobordism2)
            elif Complex.elements[k].arcs[Complex.elements[k].top +i] == Complex.elements[k].top+i+1: # source is closed but target is open, add 2 cobordisms to newRow only
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newcomps = [compcup(component) for component in cob.comps]
                    def computedeco2(decoration):
                        if decoration[magic_index+1] == 1:
                            return incrementH(decoration)
                        else:
                            return adddotonmagicindex(decoration)
                    newDecos1 = [deco for deco in cob.decos]
                    newDecos2 = [computedeco2(deco) for deco in cob.decos]
                    NewDecos1 = [deco for deco in newDecos1 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    NewDecos2 = [deco for deco in newDecos2 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    simplify_decos(NewDecos2)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], NewDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[1], bots[0], NewDecos2, newcomps)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
            else: # source and target are both open, add 1 cobordism to newRow only
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = []
                    magic_index_2 = 0 # index of component containing i+1
                    for x,comp in enumerate(cob.comps):
                        if cob.front.top+i+1 in comp:
                            magic_index_2 = x
                    if magic_index == magic_index_2: # only one component being connected via cup
                        comp = cob.comps[magic_index]
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
                        comp1 = cob.comps[magic_index]
                        comp2 = cob.comps[magic_index_2]
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
                                return incrementH(decoration)[:magic_index_2+1] + incrementH(decoration)[magic_index_2+2:]
                            elif decoration[magic_index+1] == 1 or decoration[magic_index_2+1]==1:
                                return adddotonmagicindex(decoration)[:magic_index_2+1] + adddotonmagicindex(decoration)[magic_index_2+2:]
                            else:
                                return decoration[:magic_index_2+1] + decoration[magic_index_2+2:]
                        newDecos1 = [computedeco2comps(deco) for deco in cob.decos]
                    NewDecos1 = [deco for deco in newDecos1 if deco[find_first_index(newcomps,contains_0)+1] == 0]
                    simplify_decos(NewDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], NewDecos1, newcomps)
                    newRow.append(Cobordism1)
        NewMorphisms.append(newRow)
        if Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # target is closed
            NewMorphisms.append(nextRow)
    return ChainComplex(newElements, NewMorphisms)
  
def AddPosCrossing(Complex, i):
    CapCup = AddCap(AddCup(Complex, i), i, "true")
    sourceelements = Complex.elements
    targetelements = CapCup.elements
    NewElements = sourceelements + targetelements
    TopLeft = Complex.morphisms
    BottomRight = CapCup.morphisms
    length1 = len(Complex.elements)
    length2 = len(CapCup.elements)
    TopRight = np.full((length1, length2), ZeroCob, Cobordism)
    BottomLeft = []
    for x, targetclt in enumerate(Complex.elements):
        if targetclt.arcs[targetclt.top +i] == targetclt.top+i+1: #targetclt is closed
            newRow = []
            nextRow = []
            for y, sourceclt in enumerate(Complex.elements):
                if x == y:
                    newTarget1 = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
                    newTarget2 = AddCapToCLT(AddCupToCLT(targetclt, i)[1], i, "true")
                    newcomps = components(sourceclt, newTarget1)
                    magic_index = 0
                    for x,comp in enumerate(newcomps):
                        if sourceclt.top +i in comp:
                            magic_index = x
                            break
                    decos1 = [[0] + [0 for comp in newcomps[:magic_index]] + [1] + [0 for comp in newcomps[magic_index+1:]] + [1], [1] + [0 for comp in newcomps] + [-1]] #Compute new decos via neckcutting
                    decos2 = [[0] + [0 for comp in newcomps] + [1]]
                    NewCobordism1 = Cobordism(sourceclt, newTarget1, decos1)
                    NewCobordism2 = Cobordism(sourceclt, newTarget2, decos2)
                    newRow.append(NewCobordism1)
                    nextRow.append(NewCobordism2)
                else:
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
            BottomLeft.append(newRow)
            BottomLeft.append(nextRow)
        else:
            newRow = []
            for y, sourceclt in enumerate(Complex.elements):
                if x == y:
                    newTarget = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
                    decos = [[0] + [0 for comp in components(sourceclt, newTarget)] + [1]]
                    NewCobordism = Cobordism(sourceclt, newTarget, decos)
                    newRow.append(NewCobordism)
                else:
                    newRow.append(ZeroCob)
            BottomLeft.append(newRow)
    NewMorphisms = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewComplex = ChainComplex(NewElements, NewMorphisms)
    return NewComplex

def grshiftclt(clt):
    return CLT(clt.top, clt.bot, clt.arcs, clt.pgr+1, clt.qgr+1, clt.dgr-0.5)

def grshiftcob(cob):
    newDecos = [deco[:-1] + [deco[-1]*-1] for deco in cob.decos]
    return Cobordism(grshiftclt(cob.front), grshiftclt(cob.back), newDecos, cob.comps)
    
def AddNegCrossing(Complex, i):
    complexgradings = [ [clt.pgr, clt.qgr] for clt in Complex.elements]
    targetelements = [grshiftclt(clt) for clt in Complex.elements]
    newgradings = [ [clt.pgr, clt.qgr] for clt in targetelements]
    CapCup = AddCap(AddCup(Complex, i), i)
    sourceelements = CapCup.elements
    NewElements = sourceelements + targetelements
    TopLeft = CapCup.morphisms
    BottomRight = [[grshiftcob(cob) for cob in row] for row in Complex.morphisms]
    length1 = len(sourceelements)
    length2 = len(targetelements)
    TopRight = np.full((length1, length2), ZeroCob, Cobordism)
    BottomLeft = [] 
    for x, targetclt in enumerate(Complex.elements):
        newRow = []
        for y, sourceclt in enumerate(Complex.elements):
            if sourceclt.arcs[sourceclt.top +i] == sourceclt.top+i+1: #sourceclt is closed
                if x == y:
                    newSource1 = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
                    newSource2 = AddCapToCLT(AddCupToCLT(sourceclt, i)[1], i)
                    newTarget = targetelements[x]
                    newcomps = components(newSource2, newTarget)
                    magic_index = 0
                    for z,comp in enumerate(newcomps):
                        if sourceclt.top +i in comp:
                            magic_index = z
                            break
                    decos1 = [[0] + [0 for comp in newcomps] + [1]]
                    decos2 = [[0] + [0 for comp in newcomps[:magic_index]] + [1] + [0 for comp in newcomps[magic_index+1:]] + [1]]
                    newCobordism1 = Cobordism(newSource1, newTarget, decos1) #fill out decos and comps
                    newCobordism2 = Cobordism(newSource2, newTarget, decos2)
                    newRow.append(newCobordism1)
                    newRow.append(newCobordism2)
                else:
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
            else: #sourceclt is open
                if x == y:
                    newSource = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
                    newTarget = targetelements[x]
                    decos = [[0] + [0 for comp in components(newSource, newTarget)] + [1]]
                    newCobordism = Cobordism(newSource, newTarget, decos)
                    newRow.append(newCobordism)
                else:
                    newRow.append(ZeroCob)
        BottomLeft.append(newRow)
    NewMorphisms = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewComplex = ChainComplex(NewElements, NewMorphisms)
    return NewComplex

def BNbracket(string,pos=0,neg=0,start=1):
    """compute the Bar-Natan bracket for tangle specified by 'string', which is a concatenation of words <type>+<index>, separated by '.' read from right to left, for each elementary tangle slice, read from top to bottom, where:
    <type> is equal to:
        'pos': positive crossing
        'neg': negative crossing
        'cup': cup
        'cap': cap
    <index> is the index at which the crossing, cap or cup sits. 
    'pos' and 'neg' are the numbers of positive and negative crossings. 
    The optional paramter 'start' is an integer which specifies the number of tangle ends at the top.
    E.g. 'BNbracket('cup0pos0',2)' is a (2,0)-tangle which is decomposed as a positive crossing followed by a cap.
    """
    stringlist=[[word[0:3],int(word[3:])] for word in string.split('.')]
    stringlist.reverse()
    cx=ChainComplex([CLT(start,start,[start+i for i in range(start)]+[i for i in range(start)], 0,0,0)], [[ZeroCob]])
    
    print("Computing the Bar-Natan bracket for the tangle\n\n"+string+"\n\n"+"with "+str(start)+" ends at the top, "+str(pos)+" positive crossings and "+str(neg)+" negative crossings.")
    
    for i,word in enumerate(stringlist):
        print("slice "+str(i)+": adding "+word[0]+" at index "+str(word[1])+" to tangle.", end='\r')# monitor
        #time.sleep(0.1)
        
        if word[0]=="pos":
            cx=AddPosCrossing(cx, word[1])
            cx.eliminateAll()
        
        if word[0]=="neg":
            cx=AddNegCrossing(cx, word[1])
            cx.eliminateAll()
        
        if word[0]=="cup":
            cx=AddCup(cx, word[1])
            cx.eliminateAll()
        
        if word[0]=="cap":
            cx=AddCap(cx, word[1])
    
    cx.shift_qhd(pos-2*neg,-neg,0.5*neg)
    
    print("Completed the computation successfully.       ")
    return cx


