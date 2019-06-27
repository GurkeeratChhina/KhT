import numpy as np
import math
from KhT import *
from Tangles import *
from Cobordisms import *
from Drawing import *

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

def AddCapToCLT(clt, i, grshift = "false"):
    def incrementby2(j): #increment TEI by 2 if greater than where cap is to be inserted
        if j >= clt.top+i:
            return j+2
        else:
            return j
    newarcs = [incrementby2(x) for x in clt.arcs[:clt.top+i]] + [clt.top+i +1, clt.top+i] + [incrementby2(x) for x in clt.arcs[clt.top+i:]]
    if grshift == "true":
        newgr = [j+1 for j in clt.gr]
    else:
        newgr = clt.gr
    return CLT(clt.top, clt.bot+2, newarcs, newgr)

def AddCap(Complex, i, grshift = "false"):
    """ Adds a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot """
    NewElements = [AddCapToCLT(clt, i, grshift) for clt in Complex.elements]
    length = len(Complex.elements)
    NewMorphisms = []
    for j ,row in enumerate(Complex.morphisms):
        NewRow=[]
        for k, cob in enumerate(row):
            if cob.decos == []:
                NewRow.append(ZeroCob)
            else:
                def incrementby2(j):
                    if  j >= cob.front.top+i: #increment TEI by 2 if greater than where cap is to be inserted
                        return j+2
                    else:
                        return j
                newcomps=[[incrementby2(x) for x in comp] for comp in cob.comps] + [[cob.front.top +i, cob.front.top+i+1]]         
                NewCob = Cobordism(NewElements[j], NewElements[k], [ NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos], newcomps)
                NewRow.append(NewCob)
        NewMorphisms.append(NewRow)
    return ChainComplex(NewElements, NewMorphisms)

def AddCupToCLT(clt, i): #TODO: check usages of .copy()
    """ Here 0 <= i <= clt.bot
        Returns a list of 2 clt if there is a closed component
        otherwise returns a list of 1 clt """
    newElements = []
    def decrementby2(j):#decrement TEI by 2 if greater than where cup is to be inserted
            if j >= clt.top +i:
                return j-2
            else:
                return j
    if clt.arcs[clt.top +i] == clt.top+i+1: # adding the cup makes a closed component
        newarcs = [decrementby2(x) for x in clt.arcs if x != clt.top+i and x!= clt.top+i+1]
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, [clt.gr[0]+1, clt.gr[1]]))
        newElements.append(CLT(clt.top, clt.bot -2, newarcs, [clt.gr[0]-1, clt.gr[1]]))
    else: # adding the cup doesnt make closed components
        leftend = clt.arcs[clt.top + i]
        rightend = clt.arcs[clt.top+i+1]
        newarcs = clt.arcs.copy()
        newarcs[leftend] = rightend # connecting the two arcs via the cup
        newarcs[rightend] = leftend
        newarcs1 = [decrementby2(entry) for x, entry in enumerate(newarcs) if x != clt.top+i and x != clt.top+i+1]
        newElements.append(CLT(clt.top, clt.bot -2, newarcs1, clt.gr))
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
            for x,comp in enumerate(cob.comps): #TODO: find magic index without iterating twice, dont need defining function
                if cob.front.top +i in comp:
                    magic_index = x
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
    TopLeft = Complex.morphisms
    BottomRight = CapCup.morphisms
    length1 = len(Complex.elements)
    length2 = len(CapCup.elements)
    TopRight = np.full((length1, length2), ZeroCob, Cobordism)
    BottomLeft = np. full((length2, length1), ZeroCob, Cobordism) #compute this manually
    NewMorphisms = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewElements = Complex.elements + CapCup.elements
    NewComplex = ChainComplex(NewElements, NewMorphisms)
    return NewComplex

def grshiftclt(clt):
    return CLT(clt.top, clt.bot, clt.arcs, [j+1 for j in clt.gr])

def grshiftcob(cob):
    return Cobordism(grshiftclt(cob.front), grshiftclt(cob.back), cob.decos, cob.comps)
    
def AddNegCrossing(Complex, i):
    CapCup = AddCap(AddCup(Complex, i), i)
    targetelements = [grshiftclt(clt) for clt in Complex.elements]
    sourceelements = CapCup.elements
    print("sourceelements", sourceelements)
    print("targetelements", targetelements)
    NewElements = sourceelements + targetelements
    TopLeft = CapCup.morphisms
    BottomRight = [[grshiftcob(cob) for cob in row] for row in Complex.morphisms]
    length1 = len(sourceelements)
    length2 = len(targetelements)
    TopRight = np.full((length1, length2), ZeroCob, Cobordism)
    BottomLeft = [] #np. full((length2, length1), ZeroCob, Cobordism) #compute this manually
    for x, targetclt in enumerate(Complex.elements): #maybe dont need enumerate here
        newRow = []
        for sourceclt in Complex.elements:
            fronts = iter(sourceelements)
            if sourceclt.arcs[sourceclt.top +i] == sourceclt.top+i+1: #sourceclt is closed
                print("closed")
                decos1 = []
                decos2 = []
                newCobordism1 = Cobordism(next(fronts), targetelements[x], decos1) #fill out decos and comps
                newCobordism2 = Cobordism(next(fronts), targetelements[x], decos2)
                newRow.append(newCobordism1)
                newRow.append(newCobordism2)
            else: #sourceclt is open
                print("open")
                decos = []
                newCobordism = Cobordism(next(fronts), targetelements[x], decos)
                newRow.append(newCobordism)
        BottomLeft.append(newRow)
    print(length1, length2, BottomLeft)
    NewMorphisms = np.concatenate((np.concatenate((TopLeft, TopRight), axis = 1), \
                                   np.concatenate((BottomLeft, BottomRight), axis = 1)), axis = 0)
    NewComplex = ChainComplex(NewElements, NewMorphisms)
    return NewComplex