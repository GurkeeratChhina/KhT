# -*- coding: utf-8 -*-

import math
from KhT import *
from Tangles import *
from Cobordisms import *
from Complex import *
from Drawing import *

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

# Elementary tangles and cobordisms
b=CLT(2,2,[1,0,3,2],[0,0])
c=CLT(2,2,[2,3,0,1],[0,0])
Sbc=Cobordism(b,c,[[0,0,1]])
Scb=Cobordism(c,b,[[0,0,1]])
MinusScb = Cobordism(c,b,[[0,0,-1]])
RightDc = Cobordism(c,c,[[0,0,1,1]])


def TestSet4():
    drawclt(b,"b") 
    drawclt(c,"c")
    drawcob(Sbc,"Sbc")
    drawcob(Scb,"Scb")
    drawcob((Sbc*Scb),"SSbb")
    drawcob((Scb*Sbc),"SScc")
    drawcob(RightDc, "RightDc")
    T3=CLT(4,4,[4,7,3,2,0,6,5,1],[0,0])
    T4=CLT(4,4,[6,7,3,2,5,4,0,1],[0,0])
    T5=CLT(4,4,[3,2,1,0,7,6,5,4],[0,0])
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
    
# testcobordism=Cobordism(b,c,[[2,1,2],[0,0,4]])
# drawclt(cup(10,3)*cap(10,2)*cup(10,4)*cup(8,6)+cap(2,1),"test1")
# drawclt(cap(5,5),"test2")
# drawcob(testcobordism,"testcobordism")

T1=CLT(2,4,[2,5,0,4,3,1],[0,0])
T2=CLT(2,4,[2,3,0,1,5,4],[0,0])
cob1=Cobordism(T1,T2,[[4,1,0,1]])
cob2=Cobordism(T1,T2,[[4,1,0,-3],[2,0,1,1],[1,1,1,19]])

def TestSet5():
    complex1 = ChainComplex([b,c], [[ZeroCob, ZeroCob], [Sbc ,ZeroCob]])
    complex1.ValidMorphism()

    CobRightDotMinusLeftDotVertical = Cobordism(c,c, [[0,0,1,1],[0,1,0,-1]])
    drawcob(CobRightDotMinusLeftDotVertical, "temporary1")

    complex2 = ChainComplex([b,c,c], [[ZeroCob, ZeroCob, ZeroCob],[Sbc, ZeroCob, ZeroCob],[ZeroCob, CobRightDotMinusLeftDotVertical, ZeroCob]])
    complex2.ValidMorphism()
    drawcob(Sbc*CobRightDotMinusLeftDotVertical, "temporary2")
    drawcob(Sbc*RightDc, "DottedSaddle")
   
    complex1cap0 = AddCap(complex1, 0)
    #PrettyPrintComplex(complex1cap0)
    #PrettyPrintComplex(complex1cap0,"long")
    
    DrawFourEndedChainComplex(complex2, "complex2.png")
    DrawFourEndedChainComplex(complex1, "complex1.png")
    
    complex1again = AddCup(complex1cap0, 1)
    print("printing components")
    for row in complex1again.morphisms:
        for cob in row:
            print(cob.comps)
    DrawFourEndedChainComplex(complex1again, "complex1again.png")

    RightDcMinusH = Cobordism(c,c,[[0,0,1,1], [1,0,0,-1]])

    complex3 = ChainComplex([b,c,c,c], [[ZeroCob, ZeroCob, ZeroCob, ZeroCob],[Sbc,ZeroCob, ZeroCob, ZeroCob],[ZeroCob, RightDc, ZeroCob, ZeroCob],[ZeroCob, ZeroCob, RightDcMinusH, ZeroCob]])
    DrawFourEndedChainComplex(complex3, "complex3.png")
    complex3.ValidMorphism()

    complex4 = ChainComplex([c,b,c,c], [[ZeroCob, MinusScb, RightDc, ZeroCob], [ZeroCob, ZeroCob, ZeroCob, Sbc], [ZeroCob,ZeroCob, ZeroCob, Cobordism(c,c,[[0,0,0,1]])], [ZeroCob, ZeroCob, ZeroCob, ZeroCob]])
    DrawFourEndedChainComplex(complex4, "complex4.png")
    
    

complex5 = ChainComplex([c,b,b,c,b,b], [[ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [Scb, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [ZeroCob, Cobordism(b,b, [[0,0,1,1]]), ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [Cobordism(c, c, [[1,0,0,1]]), ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [ZeroCob, Cobordism(b,b,[[1, 0,0, -1]]), ZeroCob, Scb, ZeroCob, ZeroCob],\
                                        [ZeroCob, ZeroCob, Cobordism(b,b,[[1, 0,0,1]]), ZeroCob, Cobordism(b,b,[[0,0,1,1]]), ZeroCob]])
DrawFourEndedChainComplex(complex5, "complex5.png")
complex5.ValidMorphism()

NewT1 = AddCapToCLT(T1, 1)
drawclt(NewT1, "NewT1")

NewT2 = AddCapToCLT(T2, 1)
drawclt(NewT2, "NewT2")

def newCob(cob, i):
    def incrementindex(entry):
        if  entry >= cob.front.top+i:
            return entry+2
        else:
            return entry
    newcomps=[[incrementindex(entry) for entry in comp] for comp in cob.comps]
    newcomps.append([cob.front.top +i, cob.front.top+i+1])          
    return Cobordism(AddCapToCLT(cob.front, i), AddCapToCLT(cob.back, i), [ NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos], newcomps)

drawcob(newCob(cob2, 1), "NewCob")

def TestSet0():
    BasicComplex = ChainComplex([CLT(1,1, [1,0], [0,0])], [[ZeroCob]])
    BasicCap = AddCap(BasicComplex, 0)

    drawclt(BasicCap.elements[0], "basiccap")
    Double = AddCup(BasicCap, 0)
    drawclt(Double.elements[0], "double")

    PrintComplexMorphismIntMatrix(Double)

def TestSet1():
    tempcob1 = Cobordism(CLT(1,3, [1, 0, 3, 2], [0,0]), CLT(1,3, [1, 0, 3, 2], [0,1]), [[0, 0, 1, 1]])
    tempcomplex = ChainComplex([CLT(1,3, [1, 0, 3, 2], [0,0]), CLT(1,3, [1, 0, 3, 2], [0,1])], [[ZeroCob, ZeroCob], [tempcob1, ZeroCob]])
    tempcomplexwithcup = AddCup(tempcomplex, 1)
    PrintComplexMorphismIntMatrix(tempcomplexwithcup)
    print("Decorations at 3, 0", tempcomplexwithcup.morphisms[3][0].decos)
    print("Decorations at 3, 1", tempcomplexwithcup.morphisms[3][1].decos)
    print("Grading of element at 0", tempcomplexwithcup.elements[0].gr)
    print("Grading of element at 1", tempcomplexwithcup.elements[1].gr)
    print("Grading of element at 2", tempcomplexwithcup.elements[2].gr)
    print("Grading of element at 3", tempcomplexwithcup.elements[3].gr)

def TestSet2():
    complex5cap = AddCap(complex5, 1)
    PrintComplexMorphismIntMatrix(complex5cap)
    PrintComplexMorphismDecoCompMatrix(complex5cap)
    complex5again = AddCup(complex5cap, 0)
    DrawFourEndedChainComplex(complex5again, "complex5again.png")

def TestSet3():
    CLT01 = CLT(1,3,[1,0,3,2],[0,0])
    CLT03 = CLT(1,3,[3,2,1,0],[0,1])
    BasicSaddle = Cobordism(CLT01, CLT03, [[0,0,1]])
    BasicSaddleComplex = ChainComplex([CLT01, CLT03], [[ZeroCob, ZeroCob], [BasicSaddle, ZeroCob]])
    BasicSaddleCup = AddCup(BasicSaddleComplex, 1)
    PrintComplexMorphismIntMatrix(BasicSaddleCup)
    PrintComplexMorphismDecoCompMatrix(BasicSaddleCup)
    PrettyPrintComplex(BasicSaddleCup)
    PrettyPrintComplex(BasicSaddleCup,"long")

# TestSet0()
# TestSet1()
# TestSet2()
# TestSet3()
# TestSet4()
TestSet5()


# proof of concept for matrix multiplication for matrices with customized algebra addition and multiplication\n",

class testalg(object):
	def __init__(self,x):
		self.x = x
	def __mul__(self, other):
		return testalg(2*self.x*other.x)
	def __add__(self, other):
		return testalg(2+self.x+other.x)
        

def TestAlgTest():
    a=testalg(1)
    b=testalg(2)
    c=testalg(3)
    d=testalg(4)
    A=[[a,b],[c,d]]
    print((np.tensordot(A,A, axes=(-1,-2)))[0,0].x)        
    B=[[1,2],[3,4]]
    C=[[1,1],[0,0]]
    print(np.tensordot(C,B, axes=(-1,-2)))

# TestAlgTest()

## comparing efficiency of two functions
#import timeit
#print(timeit.timeit('simplify_decos_old([[1,0,1,1],[0,0,1,2],[0,0,1,3],[0,1,1,5],[1,0,1,12]])',setup="from __main__ import simplify_decos_old", number=1000)) 
#print(timeit.timeit('simplify_decos([[1,0,1,1],[0,0,1,2],[0,0,1,3],[0,1,1,5],[1,0,1,12]])',setup="from __main__ import simplify_decos", number=1000))
## old vs new: 0.14744883000093978 vs 0.009813104006752837
print("File executed successfully")

