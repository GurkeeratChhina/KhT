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
        # return True
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

def AddCupToCLT(clt, i): #check usages of .copy()
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

def AddCup(Complex, i): # TODO: Reduce/eliminate usages of .copy(), reduce decorations
    """ Here 0 <= i <= tangle.bot -2"""
    newElements = []
    for clt in Complex.elements:
        newElements.extend(AddCupToCLT(clt, i))
    NewMorphisms = []
    for j, row in enumerate(Complex.morphisms): # j is the target clt
        newRow = []
        nextRow = []
        for k, cob in enumerate(row): # k is the source clt
            def shift(j): #shifts TEI by -2 if it is greater than where the cup i added
                if j > cob.front.top+i+1:
                    return j-2
                else:
                    return j
            magic_index = 0
            for x,comp in enumerate(cob.comps): #TODO: find magic index without iterating twice, dont need defining function
                if cob.front.top +i in comp:
                    magic_index = x
            def compcup(component): #computes the new component after adding a cup, for components containing at least 3 ends
                return [shift(j) for j in component if j != cob.front.top+i and j!= cob.front.top+i+1]
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
                    newcomps = [[shift(j) for j in comp] for comp in cob.comps if cob.front.top+i not in comp ]
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
                    simplify_decos(newDecos4)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[1], bots[0], [], newcomps)
                    Cobordism3 = Cobordism(tops[0], bots[1], newDecos3, newcomps)
                    Cobordism4 = Cobordism(tops[1], bots[1], newDecos4, newcomps)
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
                    simplify_decos(newDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[0], bots[1], newDecos2, newcomps)
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
                    simplify_decos(newDecos2)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, newcomps)
                    Cobordism2 = Cobordism(tops[1], bots[0], newDecos2, newcomps)
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
                        comp1 = comp[:min(x1, x2)] + comp[max(x1, x_2)+1:]
                        comp2 = comp[min(x1, x_2) +1:max(x1, x_2)]
                        newcomps = [[shift(j) for j in comp] for comp in cob.comps[:magic_index] +[comp1, comp2] + cob.comps[magic_index+1:]]
                        def computedeco1comps(decoration):
                            if decoration[magic_index+1] == 1:
                                return [decoration[:magic_index+1] +[1] + decoration[magic_index+1:]]
                            else:
                                return [decoration[:magic_index+1] +[1] + decoration[magic_index+1:] , \
                                        decoration[:magic_index+2] +[1] + decoration[magic_index+2:], \
                                        decoration[0]+1 + decoration[1:magic_index+1] +[0] + decoration[magic_index+1:-1] + decoration[-1]*-1]
                        for deco in cob.decos:
                            newdecos1.extend(computedeco1comps(deco))
                    else: # two seperate components being connected via cup
                        comp1 = cob.comps[magic_index]
                        comp2 = cob.comps[magic_index_2]
                        location1 = comp1.index(cob.front.top+i)
                        location2 = comp2.index(cob.front.top+i+1)
                        comp2 = (comp2[location2+1:] + comp2[:location2]) #rotates list until i+1 is at front, and removes it
                        if location1 %2 == location2%2: # if top/bot dont line up, flip comp2
                            comp2.reverse()
                        comp1 = comp1[:location1] + comp2 + comp1[location1+1:] # insert comp2 into comp1, and dont include element i
                        newcomps = [[shift(j) for j in comp] for comp in cob.comps[:magic_index] + [comp1] + cob.comps[magic_index+1:]]
                        del newcomps[magic_index_2]
                        def computedeco2comps(decoration):
                            if decoration[magic_index+1] == 1 and decoration[magic_index_2+1]==1:
                                return incrementH(decoration)[:magic_index_2+1] + incrementH(decoration)[magic_index_2+2:]
                            elif decoration[magic_index+1] == 1 or decoration[magic_index_2+1]==1:
                                return adddotonmagicindex(decoration)[:magic_index_2+1] + adddotonmagicindex(decoration)[magic_index_2+2:]
                            else:
                                return decoration[:magic_index_2+1] + decoration[magic_index_2+2:]
                        newDecos1 = [computedeco2comps(deco) for deco in cob.decos]
                    simplify_decos(newDecos1)
                    tops = AddCupToCLT(cob.front, i)
                    bots = AddCupToCLT(cob.back, i)
                    Cobordism1 = Cobordism(tops[0], bots[0], newDecos1, newcomps)
                    newRow.append(Cobordism1)
        NewMorphisms.append(newRow)
        if Complex.elements[j].arcs[clt.top +i] == clt.top+i+1: # target is closed
            NewMorphisms.append(nextRow)
    return ChainComplex(newElements, NewMorphisms)
    
def AddPosCrossing(Complex, i):
    return 0

def AddNegCrossing(Complex, i):
    return 0