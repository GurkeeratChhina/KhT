from itertools import product
import itertools as itertools
import numpy as np
import math
from KhT import *
from Tangles import *
from Cobordisms import *

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
        
        squared = flatten(np.tensordot(self.morphisms,self.morphisms, axes=(-2,-1)))
        for i in squared:
            if i.decos != []:
                raise Exception('Differential does not square to 0')

def AddCapToCLT(clt, i):
    newarcs = clt.arcs.copy()
    for j, x in enumerate(newarcs):
        if x >= clt.top + i:
            newarcs[j] += 2
    newarcs.insert(clt.top+i, clt.top+i)
    newarcs.insert(clt.top+i, clt.top+i+1)
    return CLT(clt.top, clt.bot+2, newarcs, clt.gr)

def AddCap(Complex, i):
    """ Adds a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot """
    NewElements = []
    for clt in Complex.elements:
        NewElements.append(AddCapToCLT(clt, i))
    length = len(Complex.elements)
    NewMorphisms = []
    for j ,row in enumerate(Complex.morphisms):
        NewRow=[]
        for k, cob in enumerate(row):
            if cob.decos == []:
                NewRow.append(ZeroCob)
            else:
                def incrementindex(entry):
                    if  entry >= cob.front.top+i:
                        return entry+2
                    else:
                        return entry
                newcomps=[[incrementindex(entry) for entry in comp] for comp in cob.comps]
                newcomps.append([cob.front.top +i, cob.front.top+i+1])          
                NewCob = Cobordism(NewElements[j], NewElements[k], [ NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos], newcomps)
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

def AddCup(Complex, i): # WIP
    """ Here 0 <= i <= tangle.bot -2"""
    newElements = []
    for clt in Complex.elements:
        newElements.extend(AddCupToCLT(clt, i))
    NewMorphisms = []
    for j, row in enumerate(Complex.morphisms): # j is the target clt
        newRow = []
        nextRow = []
        for k, cob in enumerate(row): # k is the source clt
            if Complex.elements[k].arcs[Complex.elements[k].top +i] == Complex.elements[k].top+i+1 \
                and Complex.elements[j].arcs[Complex.elements[j].top +i] == Complex.elements[j].top+i+1: # source and target are both closed, add 2 cobordisms to newRow and nextRow
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = []
                    newDecos2 = []
                    newDecos3 = []
                    newDecos4 = []
                    magic_index = cob.comps.index([cob.front.top+i, cob.front.top+i+1]) #computes index of closed component
                    for deco in cob.decos:
                        decocopy = deco.copy()
                        del decocopy[magic_index + 1] # removing closed component
                        if deco[magic_index + 1] == 0: # no dot on component already prior
                            newDecos1.append(decocopy)
                            newDecos4.append(decocopy)
                        else: # is dot on component already
                            newDecos3.append(decocopy)
                            Hcopy = decocopy.copy()
                            Hcopy[0] += 1
                            newDecos4.append(Hcopy)
                    simplify_decos(newDecos4)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1)
                    Cobordism2 = Cobordism(tops[1], bots[0], newDecos2)
                    Cobordism3 = Cobordism(tops[0], bots[1], newDecos3)
                    Cobordism4 = Cobordism(tops[1], bots[1], newDecos4)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
                    nextRow.append(Cobordism3)
                    nextRow.append(Cobordism4)
            elif Complex.elements[j].arcs[Complex.elements[j].top +i] == Complex.elements[j].top+i+1: # source is open but target is closed, add 1 cobordism to each
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    nextRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = []
                    newDecos2 = []
                    magic_index = 0
                    compscopy = cob.comps.copy()
                    for x,comp in enumerate(compscopy):
                        if cob.front.top +i in comp:
                            magic_index = x
                            compscopy[x].remove(cob.front.top+i)
                            compscopy[x].remove(cob.front.top+i+1)
                        for y, end in enumerate(comp):
                            if end > cob.front.top+i+1:
                                comp[y] -= 2
                    for deco in cob.decos:
                        decocopy = deco.copy()
                        if deco[magic_index+1] == 1: # dot on component already
                            newDecos2.append(decocopy)
                        else: # dot not on component already
                            newDecos2.append(decocopy)
                            dotcopy = decocopy.copy()
                            dotcopy[magic_index+1] = 1
                            newDecos1.append(dotcopy)
                            Hcopy = decocopy.copy()
                            Hcopy[0] += 1
                            Hcopy[-1] *= -1
                            newDecos1.append(Hcopy)
                    simplify_decos(newDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, compscopy)
                    Cobordism2 = Cobordism(tops[0], bots[1], newDecos2, compscopy)
                    newRow.append(Cobordism1)
                    nextRow.append(Cobordism2)
            elif Complex.elements[k].arcs[Complex.elements[k].top +i] == Complex.elements[k].top+i+1: # source is closed but target is open, add 2 cobordisms to newRow only
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = []
                    newDecos2 = []
                    compscopy = cob.comps.copy()
                    magic_index = 0
                    for x,comp in enumerate(compscopy):
                        if cob.front.top +i in comp:
                            magic_index = x
                            compscopy[x].remove(cob.front.top+i)
                            compscopy[x].remove(cob.front.top+i+1)
                        for y, end in enumerate(comp):
                            if end > cob.front.top+i+1:
                                comp[y] -= 2
                    for deco in cob.decos:
                        decocopy = deco.copy()
                        if deco[magic_index+1] == 1: # dot on component already
                            newDecos1.append(decocopy)
                            dotcopy = decocopy.copy()
                            Hcopy = decocopy.copy()
                            Hcopy[0] += 1
                            newDecos2.append(Hcopy)
                        else: # dot not on component already
                            newDecos1.append(decocopy)
                            dotcopy = decocopy.copy()
                            dotcopy[magic_index+1] = 1
                            newDecos2.append(dotcopy)
                    simplify_decos(newDecos2)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, compscopy)
                    Cobordism2 = Cobordism(tops[1], bots[0], newDecos2, compscopy)
                    newRow.append(Cobordism1)
                    newRow.append(Cobordism2)
            else: # source and target are both open, add 1 cobordism to newRow only
                if cob.decos == []: # is the 0 cob
                    newRow.append(ZeroCob)
                else: # is not the 0 cob
                    newDecos1 = [] #TODO: Fill out decos
                    magic_index_1 = 0 # index of component containing i
                    magic_index_2 = 0 # index of component containing i+1
                    compscopy = cob.comps.copy()
                    for x,comp in enumerate(compscopy):
                        if cob.front.top +i in comp:
                            magic_index_1 = x
                        if cob.front.top+i+1 in comp:
                            magic_index_2 = x
                    if magic_index_1 == magic_index_2: # only one component being connected via cup
                        comp = compscopy[magic_index_1].copy()
                        x1 = comp.index(cob.front.top +i)
                        x2 = comp.index(cob.front.top +i + 1)
                        comp1 = comp[:min(x1, x2)] + comp[max(x1, x_2)+1:]
                        comp2 = comp[min(x1, x_2) +1:max(x1, x_2)]
                        compscopy = compscopy[:magic_index_1] +[comp1, comp2] + compscopy[magic_index_1+1:]
                        for x, comp in enumerate(compscopy):
                            for y, end in enumerate(comp):
                                if  end > cob.front.top+i+1:
                                    comp[y] -= 2
                        for deco in cob.decos:
                            decocopy = deco.copy()
                            if decocopy[magic_index_1+1] == 1: # dot on the component we want to cut the neck of
                                decocopy.insert(magic_index_1 +1, 1)
                                newDecos1.append(decocopy)
                            else: # no dot on component we want to cut the neck of
                                decodotleft = decocopy.copy()
                                decodotright = decocopy.copy()
                                decominusH = decocopy.copy()
                                decodotleft.insert(magic_index_1 +1, 1)
                                decodotleft.insert(magic_index_1 +2, 1)
                                decominusH.insert(magic_index_1 +1, 0)
                                decominusH[0] += 1
                                decominusH[-1] *= -1
                                newDecos1.append(decodotleft)
                                newDecos1.append(decodotright)
                                newDecos1.append(decominusH)
                            
                    else: # two seperate components being connected via cup
                        comp1 = compscopy[magic_index_1]
                        comp2 = compscopy[magic_index_2]
                        location1 = comp1.index(cob.front.top+i)
                        location2 = comp2.index(cob.front.top+i+1)
                        Alongtop1 = location1%2 #0 if the arc immediately after i+1 is along top, 1 otherwise
                        Alongtop2 = location2%2 #0 if the arc immediately after i+1 is along top, 1 otherwise
                        comp2 = (comp2[-1*location2:] + comp2[:-1*location2]) #rotates list until i+1 is at front
                        del comp2[0] #remove i+1
                        if location1 != location2: # if top/bot dont line up, flip comp2
                            comp2.reverse()
                        comp1 = comp1[:location1] + comp2 + comp1[location1+1:] # insert comp2 into comp1, and dont include element i
                        del compscopy[magic_index_2]
                        for x, comp in enumerate(compscopy):
                            for y, end in enumerate(comp):
                                if end > cob.front.top+i+1:
                                    comp[y] -= 2
                        for deco in cob.decos:
                            decocopy = deco.copy()
                            if decocopy[magic_index_1 + 1] == 1 and decocopy[magic_index_2 + 1] == 1: # dot on both components
                                decocopy[0] += 1
                            elif decocopy[magic_index_1 + 1] == 1 or decocopy[magic_index_2 + 1] == 1: # dot on exactly one of the two components
                                decocopy[magic_index_1 +1] = 1
                            del decocopy[magic_index_2 +1]
                            newDecos1.append(decocopy)
                    simplify_decos(newDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, compscopy)
                    newRow.append(Cobordism1)
        NewMorphisms.append(newRow)
        if Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # target is closed
            NewMorphisms.append(nextRow)
    return ChainComplex(newElements, NewMorphisms)
    
def AddPosCrossing(Complex, i):
    return 0

def AddNegCrossing(Complex, i):
    return 0