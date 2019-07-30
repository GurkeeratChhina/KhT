import math
from KhT import *
from Tangles import *
from Complex import *
from Drawing import *
from Cobordisms import *
from BNComplexes import *

class Tangle(object):
    __slots__ = 'slices', 'strlist'
    
    def __init__(self, string):
        self.slices = string
        
    def vertical_sum(self, other):
        return Tangle("cup2." + other.slices + "." + self.slices)    
    
    def shift(self, n):
        def shift_word(word):
            return word[0:3] + str(int(word[3:]) + n)
        splitlist = [shift_word(word) for word in self.slices.split('.')]
        newstring = ""
        for word in splitlist[:-1]:
            newstring += word + "."
        newstring += splitlist[-1]
        return Tangle(newstring)
    
    def horizontal_sum(self, other):
        return Tangle("cup1." + other.shift(2).slices + "." + self.slices)
    
    def toReduced_BNComplex(self, max_iter = 1000, start = 1, field = 2, options = "unsafe"):
        pos = 0 #TODO: Compute positive and negative crossings from orientation
        neg = 0
        stringlist=[[word[0:3],int(word[3:])] for word in self.slices.split('.')]
        stringlist.reverse()
        cx=ChainComplex([CLT(start,start,[start+i for i in range(start)]+[i for i in range(start)], 0,0,0)], [[ZeroCob]])
        print("Computing the Bar-Natan bracket for the tangle\n\n"+self.slices+"\n\n"+"with "+str(start)+" ends at the top, "+str(pos)+\
              " positive crossings and "+str(neg)+" negative crossings.")
        ends = start
        for i,word in enumerate(stringlist):
            # PrettyPrintComplex(cx, "old long")
            print("slice "+str(i)+"/"+str(len(stringlist))+": adding "+word[0]+" at index "+str(word[1])+" to tangle. ("+str(len(cx.elements))+" objects)", end='\n')# monitor \n ->\r
            #time.sleep(0.1)
            if word[0]=="pos":
                cx=AddPosCrossing(cx, word[1])
                if options=="safe": cx.ValidMorphism()
                cx.eliminateAll()
                if options=="safe": cx.ValidMorphism()
                
                
            if word[0]=="neg":
                cx=AddNegCrossing(cx, word[1])
                if options=="safe": cx.ValidMorphism()
                cx.eliminateAll()
                if options=="safe": cx.ValidMorphism()
            
            if word[0]=="cup":
                cx=AddCup(cx, word[1])
                if options=="safe": cx.ValidMorphism()
                cx.eliminateAll()
                if options=="safe": cx.ValidMorphism()
                ends -= 2
            
            if word[0]=="cap":
                cx=AddCap(cx, word[1])
                if options=="safe": cx.ValidMorphism()
                ends += 2
            
            if ends == 3:
                Draw_1_3_ChainComplex(cx, self.slices+"slice_" + str(i) + "_pre_cleanup.svg")
                print("cleaning up")
                BNcx = CobComplex2BNComplex(cx, field)
                BNcx.clean_up(max_iter)
                # cx = BNComplex2CobComplex(BNcx)
                # Draw_1_3_ChainComplex(cx, self.slices+"slice_" + str(i) + "_post_cleanup.svg")
                print("finished cleaning up")
                
            
        cx.shift_qhd(pos-2*neg,-neg,0.5*neg)
        print("Completed the computation successfully.                                              ")
        BN_complex= CobComplex2BNComplex(cx, field)
        BN_complex.clean_up(max_iter)
        return BN_complex
