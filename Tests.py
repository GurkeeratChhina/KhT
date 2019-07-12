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
from KhT import *
from Tangles import *
from Cobordisms import *
from Complex import *
from Drawing import *
from BNComplexes import *
from fractions import Fraction

# create basic crossingless tangles
def cup_alt(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cup"""
    return CLT(n+1,n-1,\
               [n+1+j for j in range(i-1)]+\
               [i,i-1]+\
               [n+j for j in range(i,n)]+\
               [j for j in range(i-1)]+\
               [j+2 for j in range(i-1,n-1)]\
               ,0, 0, 0)
def parallel(n):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cap"""
    return CLT(n,n,[n+j for j in range(n)]+[j for j in range(n)],0, 0, 0)
def cup(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cap"""
    return parallel(i-1)+CLT(2,0,[1,0],0, 0, 0)+parallel(n-i)
def cap(n,i):
    """Create a CLT with n strands of which all are parallel except for the ith which is a cup"""
    return parallel(i-1)+CLT(0,2,[1,0],0, 0, 0)+parallel(n-i)

# Elementary tangles and cobordisms
b=CLT(2,2,[1,0,3,2],0,0,0)
c=CLT(2,2,[2,3,0,1],0,0,0)
Sbc=Cobordism(b,c,[[0,0,1]])
Scb=Cobordism(c,b,[[0,0,1]])
MinusScb = Cobordism(c,b,[[0,0,-1]])
RightDc = Cobordism(c,c,[[0,0,1,1]])
   
# testcobordism=Cobordism(b,c,[[2,1,2],[0,0,4]])
# drawclt(cup(10,3)*cap(10,2)*cup(10,4)*cup(8,6)+cap(2,1),"test1")
# drawclt(cap(5,5),"test2")
# drawcob(testcobordism,"testcobordism")

T1=CLT(2,4,[2,5,0,4,3,1],0,0,0)
T2=CLT(2,4,[2,3,0,1,5,4],0,0,0)
cob1=Cobordism(T1,T2,[[4,1,0,1]])
cob2=Cobordism(T1,T2,[[4,1,0,-3],[2,0,1,1],[1,1,1,19]])

#NewT1 = AddCapToCLT(T1, 1)
#drawclt(NewT1, "NewT1")

#NewT2 = AddCapToCLT(T2, 1)
#drawclt(NewT2, "NewT2")

def newCob(cob, i):
    def incrementindex(entry):
        if  entry >= cob.front.top+i:
            return entry+2
        else:
            return entry
    newcomps=[[incrementindex(entry) for entry in comp] for comp in cob.comps]
    newcomps.append([cob.front.top +i, cob.front.top+i+1])          
    return Cobordism(AddCapToCLT(cob.front, i), AddCapToCLT(cob.back, i), [ NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos], newcomps)

#drawcob(newCob(cob2, 1), "NewCob")

def TestSet0():
    BasicComplex = ChainComplex([CLT(1,1, [1,0], 0,0,0)], [[ZeroCob]])
    BasicCap = AddCap(BasicComplex, 0)

    drawclt(BasicCap.elements[0], "basiccap")
    Double = AddCup(BasicCap, 0)
    drawclt(Double.elements[0], "double")
    # PrintComplexMorphismIntMatrix(Double)
    return 0

def TestSet1():
    tempcob1 = Cobordism(CLT(1,3, [1, 0, 3, 2], 0,0,0), CLT(1,3, [1, 0, 3, 2], 1,2,0), [[0, 0, 1, 1]])
    tempcomplex = ChainComplex([CLT(1,3, [1, 0, 3, 2], 0,0,0), CLT(1,3, [1, 0, 3, 2], 1,2,0)], [[ZeroCob, ZeroCob], [tempcob1, ZeroCob]])
    tempcomplexwithcup = AddCup(tempcomplex, 1)
    return 0

def TestSet2():
    complex5 = ChainComplex([c,b,b,c,b,b], [[ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [Scb, ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [ZeroCob, Cobordism(b,b, [[0,0,1,1]]), ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [Cobordism(c, c, [[1,0,0,1]]), ZeroCob, ZeroCob, ZeroCob, ZeroCob, ZeroCob],\
                                        [ZeroCob, Cobordism(b,b,[[1, 0,0, -1]]), ZeroCob, Scb, ZeroCob, ZeroCob],\
                                        [ZeroCob, ZeroCob, Cobordism(b,b,[[1, 0,0,1]]), ZeroCob, Cobordism(b,b,[[0,0,1,1]]), ZeroCob]])
    DrawFourEndedChainComplex(complex5, "complex5.png")
    complex5.ValidMorphism()

    complex5cap = AddCap(complex5, 1)
    # PrintComplexMorphismIntMatrix(complex5cap)
    # PrintComplexMorphismDecoCompMatrix(complex5cap)
    complex5again1 = AddCup(complex5cap, 0)
    DrawFourEndedChainComplex(complex5again1, "complex5again1.png")
    complex5again2 = AddCup(AddCap(complex5, 0), 1)
    DrawFourEndedChainComplex(complex5again2, "complex5again2.png")
    
    complex5double = AddCup(AddCap(complex5, 0),0)
    DrawFourEndedChainComplex(complex5double, "complex5double.png")
    
    complex5cup = AddCup(complex5, 0)
    print("complex5cup")
    PrettyPrintComplex(complex5cup, "old long")

def TestSet3():
    CLT01 = CLT(1,3,[1,0,3,2],0,0,0)
    CLT03 = CLT(1,3,[3,2,1,0],1,1,0.5)
    BasicSaddle = Cobordism(CLT01, CLT03, [[0,0,1]])
    BasicSaddleComplex = ChainComplex([CLT01, CLT03], [[ZeroCob, ZeroCob], [BasicSaddle, ZeroCob]])
    BasicSaddleCup = AddCup(BasicSaddleComplex, 1)
    # PrintComplexMorphismIntMatrix(BasicSaddleCup)
    # PrintComplexMorphismDecoCompMatrix(BasicSaddleCup)
    # PrettyPrintComplex(BasicSaddleCup)
    print("Basic Saddle with Cup")
    PrettyPrintComplex(BasicSaddleCup,"old long")
    return 0

def TestSet4():
    drawclt(b,"b") 
    drawclt(c,"c")
    drawcob(Sbc,"Sbc")
    drawcob(Scb,"Scb")
    drawcob((Sbc*Scb),"SSbb")
    drawcob((Scb*Sbc),"SScc")
    drawcob(RightDc, "RightDc")
    T3=CLT(4,4,[4,7,3,2,0,6,5,1],0,0,0)
    T4=CLT(4,4,[6,7,3,2,5,4,0,1],0,0,0)
    T5=CLT(4,4,[3,2,1,0,7,6,5,4],0,0,0)
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
 
def TestSet5():
    complex1 = ChainComplex([b,c], [[ZeroCob, ZeroCob], [Sbc ,ZeroCob]])
    complex1.ValidMorphism()
   
    complex1cap0 = AddCap(complex1, 0)
    print("complex1cap0")
    PrettyPrintComplex(complex1cap0,"old long")
    
    DrawFourEndedChainComplex(complex1, "complex1.png")
    
    complex1again = AddCup(complex1cap0, 1)
    DrawFourEndedChainComplex(complex1again, "complex1again.png")
    
    complex1double = AddCup(complex1cap0, 0)
    DrawFourEndedChainComplex(complex1double, "complex1double.png")
    
    AddPosCrossing(complex1, 0)
    
def TestSet6():
    RightDcMinusH = Cobordism(c,c,[[0,0,1,1], [1,0,0,-1]])
    complex3 = ChainComplex([b,c,c,c], [[ZeroCob, ZeroCob, ZeroCob, ZeroCob],[Sbc,ZeroCob, ZeroCob, ZeroCob],[ZeroCob, RightDc, ZeroCob, ZeroCob],[ZeroCob, ZeroCob, RightDcMinusH, ZeroCob]])
    DrawFourEndedChainComplex(complex3, "complex3.png")
    complex3.ValidMorphism()

    complex4 = ChainComplex([c,b,c,c], [[ZeroCob, ZeroCob, ZeroCob, ZeroCob], [MinusScb, ZeroCob, ZeroCob, ZeroCob], [RightDc,ZeroCob, ZeroCob, ZeroCob], [ZeroCob, Sbc, Cobordism(c,c,[[0,0,0,1]]), ZeroCob]])
    DrawFourEndedChainComplex(complex4, "complex4.png")
    # This should fail:
    # complex4.ValidMorphism()
    
    print("complex3")
    PrettyPrintComplex(complex3, "old long")
    complex3cap0 = AddCap(complex3, 0)
    print("complex3cap")
    PrettyPrintComplex(complex3cap0, "old long")
    complex3again1 = AddCup(complex3cap0, 1)
    complex3again2 = AddCup(AddCap(complex3, 1), 0)
    complex3again3 = AddCup(AddCap(complex3, 1), 2)
    DrawFourEndedChainComplex(complex3again1, "complex3again1.png")
    DrawFourEndedChainComplex(complex3again2, "complex3again2.png")
    DrawFourEndedChainComplex(complex3again3, "complex3again3.png")
    
    complex3double1 = AddCup(AddCap(complex3, 0), 0)
    complex3double2 = AddCup(AddCap(complex3, 1), 1)
    DrawFourEndedChainComplex(complex3double1, "complex3double1.png")
    DrawFourEndedChainComplex(complex3double2, "complex3double2.png")
    
    complex4again1 = AddCup(AddCap(complex4, 0), 1)
    complex4again2 = AddCup(AddCap(complex4, 1), 0)
    DrawFourEndedChainComplex(complex4again1, "complex4again1.png")
    DrawFourEndedChainComplex(complex4again2, "complex4again2.png")
    
    complex6 = ChainComplex([c,b], [[ZeroCob, ZeroCob],[Scb, ZeroCob]])
    complex6cup = AddCup(complex6, 0)
    print("complex6cup")
    PrettyPrintComplex(complex6cup, "old long")
    
def TestSet7():
    CobRightDotMinusLeftDotVertical = Cobordism(c,c, [[0,0,1,1],[0,1,0,-1]])
    # drawcob(CobRightDotMinusLeftDotVertical, "temporary1")

    complex2 = ChainComplex([b,c,c], [[ZeroCob, ZeroCob, ZeroCob],[Sbc, ZeroCob, ZeroCob],[ZeroCob, CobRightDotMinusLeftDotVertical, ZeroCob]])
    complex2.ValidMorphism()
    # drawcob(Sbc*CobRightDotMinusLeftDotVertical, "temporary2")
    # drawcob(Sbc*RightDc, "DottedSaddle")
    
    DrawFourEndedChainComplex(complex2, "complex2.png")
    
    complex2cap1 = AddCap(complex2, 1)
    print("complex2cap1")
    PrettyPrintComplex(complex2cap1,"old long")
    complex2again = AddCup(complex2cap1, 0)
    DrawFourEndedChainComplex(complex2again, "complex2again.png")
    complex2cap0 = AddCap(complex2, 0)
    print("complex2cap0")
    PrettyPrintComplex(complex2cap0)
    complex2again2 = AddCup(complex2cap0, 1)
    print("complex2again2")
    PrettyPrintComplex(complex2again2, "old long")
    DrawFourEndedChainComplex(complex2again2, "complex2again2.png")
    
    PrettyPrintComplex(complex2cap0, "old long")
    complex2double = AddCup(complex2cap0, 0)
    PrettyPrintComplex(complex2double, "old long")
    DrawFourEndedChainComplex(complex2double, "complex2double.png")
    
    complex2double2 = AddCup(AddCap(complex2, 1), 1)
    DrawFourEndedChainComplex(complex2double2, "complex2double2.png")
    print("complex2cap")
    PrettyPrintComplex(AddCap(complex2, 1), "old long")
    print("complex2double")
    PrettyPrintComplex(complex2double2, "old long")
    
    print("complex2")
    PrettyPrintComplex(complex2, "old long")
    complex2cup = AddCup(complex2, 0)
    print("complex2cup")
    PrettyPrintComplex(complex2cup, "old long")
    print("complex2again")
    PrettyPrintComplex(complex2, "old long")

def TestSet8():
    TangleC = CLT(2,2, [2,3,0,1], 0,0,0)
    TangleB = CLT(2,2, [1,0,3,2], 0,0,0)
    BasicComplex1 = ChainComplex([TangleC], [[ZeroCob]])
    BasicComplex2 = ChainComplex([TangleB], [[ZeroCob]])
    
    temp1 = AddCap(BasicComplex1,1)
    temp2 = AddNegCrossing(temp1,0)
    temp3 = AddNegCrossing(temp2,2)
    TwoNegCrossing = AddCup(temp3,1)
    DrawFourEndedChainComplex(TwoNegCrossing, "TwoNegCrossing.png")
    print("TwoNegCrossing")
    PrettyPrintComplex(TwoNegCrossing, "old long")
    TwoNegCrossing.eliminateAll()
    DrawFourEndedChainComplex(TwoNegCrossing, "TwoNegCrossingReduced.png")
    print("TwoNegCrossingReduced")
    PrettyPrintComplex(TwoNegCrossing, "old long")
    
    ClosedOpenNeg = AddNegCrossing(BasicComplex2, 0)
    DrawFourEndedChainComplex(ClosedOpenNeg, "ClosedOpenNeg.png")
    print("ClosedOpenNeg")
    PrettyPrintComplex(ClosedOpenNeg, "old long")
    
    temp4 = AddCap(BasicComplex1,1)
    temp5 = AddPosCrossing(temp4,0)
    temp6 = AddPosCrossing(temp5,2)
    TwoPosCrossing = AddCup(temp6,1)
    DrawFourEndedChainComplex(TwoPosCrossing, "TwoPosCrossing.png")
    print("TwoPosCrossing")
    PrettyPrintComplex(TwoPosCrossing, "old long")
    
    ClosedOpenPos = AddPosCrossing(BasicComplex2, 0)
    DrawFourEndedChainComplex(ClosedOpenPos, "ClosedOpenPos.png")
    print("ClosedOpenPos")
    PrettyPrintComplex(ClosedOpenPos, "old long")
    
def TestSet9():
    TangleC = CLT(2,2, [2,3,0,1], 0,0,0)
    BasicComplex = ChainComplex([TangleC], [[ZeroCob]])
    temp1 = AddCap(BasicComplex, 1)
    temp2 = AddNegCrossing(temp1, 0)
    temp2.eliminateAll()
    print("temp2")
    PrettyPrintComplex(temp2, "old long")
    temp3 = AddNegCrossing(temp2, 0)
    print("temp3")
    PrettyPrintComplex(temp3, "old long")
    temp3.eliminateAll()
    print("temp3eliminated")
    PrettyPrintComplex(temp3, "old long")
    temp4 = AddNegCrossing(temp3, 0)
    print("temp4")
    PrettyPrintComplex(temp4, "old long")
    temp4.eliminateAll()
    print("temp4eliminate")
    PrettyPrintComplex(temp4, "old long")
    temp5 = AddPosCrossing(temp4, 2)
    temp5.eliminateAll()
    temp6 = AddPosCrossing(temp5, 2)
    temp6.eliminateAll()
    temp7 = AddCup(temp6, 1)
    temp7.eliminateAll()
    DrawFourEndedChainComplex(temp7, "-2_3_pretzel.png")

def TestSet10():
    cx=BNbracket("cup1.pos2.pos2.neg0.neg0.neg0.cap1",2)
    DrawFourEndedChainComplex(cx, "-2_3_pretzel.png")
    
def TestSet11():
    tangle="cup1.neg2.neg2.neg2.cap3.pos0.pos0.cap1"
    drawtangle(tangle,"2m3pt_diagram_slices","slices")
    drawtangle(tangle,"2m3pt_diagram_plain","plain")
    cx=BNbracket(tangle,0,5)# 0 positive twists, 5 negative twists. 
    #cx=BNbracket("cup1.pos2.pos2.pos2.cap3.neg0.neg0.cap1") #mirror
    #PrettyPrintComplex(cx)
    BNcx=CobComplex2BNComplex(cx)
    DrawBNComplex(BNcx, "2m3pt_redBN_before_cleanupX.svg","delta_h")
    #PrettyPrintBNComplex(BNcx)
    BNcx.clean_up()
    DrawBNComplex(BNcx, "2m3pt_redBN_after_cleanupX.svg","delta_h")
    #PrettyPrintBNComplex(BNcx)

def TestSet12():
    complex_4m5pt = BNbracket("cup1.neg2.neg2.neg2.neg2.neg2.cap3.pos0.pos0.pos0.pos0.cap1")
    BN_complex_4m5pt = CobComplex2BNComplex(complex_4m5pt)
    BN_complex_4m5pt.clean_up()
    DrawBNComplex(BN_complex_4m5pt, "4m5pt_redBN_after_cleanup.svg","index_qh")
    
def TestSet13():
    n = 3
    pretzeltangle_n = "cup1."
    for j in range(n):
        pretzeltangle_n += "neg2.neg2."
    pretzeltangle_n += "neg2.cap3."
    for j in range(n):
        pretzeltangle_n += "pos0.pos0."
    pretzeltangle_n += "cap1"
    drawtangle(pretzeltangle_n,str(2*n)+"m"+str(2*n+1)+"pt_diagram_slices","slices",1)
    drawtangle(pretzeltangle_n,str(2*n)+"m"+str(2*n+1)+"pt_diagram_plain","plain",1) 
    complex_ptn = BNbracket(pretzeltangle_n)
    BN_complex_ptn = CobComplex2BNComplex(complex_ptn)
    DrawBNComplex(BN_complex_ptn,str(2*n)+"m"+str(2*n+1)+"pt_redBN_before_cleanup.svg","qh")
    BN_complex_ptn.clean_up()
    DrawBNComplex(BN_complex_ptn,str(2*n)+"m"+str(2*n+1)+"pt_redBN_after_cleanup.svg","qh")

def TestSet14():
    n=2
    pretzeltangle_n = "cup1."
    for j in range(n):
        pretzeltangle_n += "neg2.neg2."
    pretzeltangle_n += "neg2."
    for j in range(n):
        pretzeltangle_n += "pos0.pos0."
    pretzeltangle_n += "cap1" 
    drawtangle(pretzeltangle_n,str(2*n)+"m"+str(2*n+1)+"pt_diagram_slicesX","slices",2)
    drawtangle(pretzeltangle_n,str(2*n)+"m"+str(2*n+1)+"pt_diagram_plainX","plain",2)  
    complex_ptn = BNbracket(pretzeltangle_n,0,0,2)
    DrawFourEndedChainComplex(complex_ptn, "pretzeltangle_n.png")
    
def TestSet15():
    n = 0
    tangle_pos_0 = "cup1.neg2.pos0.pos0.neg1."
    for j in range(n):
        tangle_pos_0 += "pos2."
    tangle_pos_0 += "neg1.pos2.cap3.neg0.neg0.neg0.cap1"
    complex_pos_0 = BNbracket(tangle_pos_0)
    BN_complex_pos_0 = CobComplex2BNComplex(complex_pos_0)
    BN_complex_pos_0.clean_up()
    DrawBNComplex(BN_complex_pos_0, "tangle_BN_0_after_cleanup.svg", "index_qh")
    
def TestSet16():
    n = 1
    tangle_pos_1 = "cup1.neg2.pos0.pos0.neg1."
    for j in range(n):
        tangle_pos_1 += "pos2."
    tangle_pos_1 += "neg1.pos2.cap3.neg0.neg0.neg0.cap1"
    complex_pos_1 = BNbracket(tangle_pos_1)
    BN_complex_pos_1 = CobComplex2BNComplex(complex_pos_1)
    BN_complex_pos_1.clean_up(2000)
    DrawBNComplex(BN_complex_pos_1, "tangle_BN_1_after_cleanup.svg", "index_qh")

def TestSet17():
    n = 3
    tangle_pos_3 = "cup1.neg2.pos0.pos0.neg1."
    for j in range(n):
        tangle_pos_3 += "pos2."
    tangle_pos_3 += "neg1.pos2.cap3.neg0.neg0.neg0.cap1"
    complex_pos_3 = BNbracket(tangle_pos_3)
    BN_complex_pos_3 = CobComplex2BNComplex(complex_pos_3)
    BN_complex_pos_3.clean_up()
    DrawBNComplex(BN_complex_pos_3, "tangle_BN_3_after_cleanup.svg", "index_qh")

def TestSet18():
    n = 4
    tangle_pos_4 = "cup1.neg2.pos0.pos0.neg1."
    for j in range(n):
        tangle_pos_4 += "pos2."
    tangle_pos_4 += "neg1.pos2.cap3.neg0.neg0.neg0.cap1"
    complex_pos_4 = BNbracket(tangle_pos_4)
    complex_pos_4.ValidMorphism()
    BN_complex_pos_4 = CobComplex2BNComplex(complex_pos_4)
    BN_complex_pos_4.clean_up()
    DrawBNComplex(BN_complex_pos_4, "tangle_BN_4_after_cleanup.svg", "index_qh")

def TestSet19():
    n = 2
    tangle_neg_4 = "cup1.neg2.pos0.pos0.neg1."
    for j in range(n):
        tangle_neg_4 += "neg2."
    tangle_neg_4 += "neg1.pos2.cap3.neg0.neg0.neg0.cap1"
    drawtangle(tangle_neg_4,"test","slices",1)
    complex_neg_4 = BNbracket(tangle_neg_4)
    #complex_neg_4.ValidMorphism()
    BN_complex_neg_4 = CobComplex2BNComplex(complex_neg_4,2)
    BN_complex_neg_4.ValidMorphism()
    BN_complex_neg_4.clean_up(2,1000)
    DrawBNComplex(BN_complex_neg_4, "tangle_BN_neg_4_after_cleanup.svg", "index_qh")
    
# TestSet0()
# TestSet1()
# TestSet2()
# TestSet3()
# TestSet4()
# TestSet5()
# TestSet6()
# TestSet7()
# TestSet8()
# TestSet9()
# TestSet10()
# TestSet11()
# TestSet12()
# TestSet13()
# TestSet14()
# TestSet15()
# TestSet16()
# TestSet17() # Doesn't seem to work
# TestSet18() # Also doesn't seem to work
TestSet19() # Also doesn't seem to work


# TestAlgTest()

## comparing efficiency of two functions
#import timeit
#print(timeit.timeit('simplify_decos_old([[1,0,1,1],[0,0,1,2],[0,0,1,3],[0,1,1,5],[1,0,1,12]])',setup="from __main__ import simplify_decos_old", number=1000)) 
#print(timeit.timeit('simplify_decos([[1,0,1,1],[0,0,1,2],[0,0,1,3],[0,1,1,5],[1,0,1,12]])',setup="from __main__ import simplify_decos", number=1000))
## old vs new: 0.14744883000093978 vs 0.009813104006752837
print("File executed successfully")
