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
    
    def toReduced_BNComplex(self, max_iter = 100, start = 1, field = 2, options = "unsafe"):
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
                BNcx = CobComplex2BNComplex(cx, field)
                BNcx.clean_up(max_iter)
                cx = BNComplex2CobComplex(BNcx)

                
            
        cx.shift_qhd(pos-2*neg,-neg,0.5*neg)
        print("Completed the computation successfully.                                              ")
        BN_complex= CobComplex2BNComplex(cx, field)
        BN_complex.clean_up(max_iter)
        return BN_complex

def ValidCup(word, i): #TODO: check if adding cup at index i is valid
    return True
def ValidPos(word, i): #TODO: check if adding pos at index i is valid
    return True
def ValidNeg(word, i): #TODO: check if adding neg at index i is valid
    return True


def GenerateTangleWords(max_length):
    remaining_slices = max_length-1
    currentListofWords = [ ["cap1", -1] ]
    while remaining_slices > 1:
        nextListofWords = []
        for word in currentListofWords:
            ends = 5 + 2*word[1]
            if abs(word[1]) > remaining_slices -1:
                continue #Do nothing, the tangle won't be a 4-ended tangle no matter what we do                
            elif abs(word[1]) == remaining_slices -1 and word[1] < 0: 
                continue #Do nothing as can only add caps to get 4 ends in time, and that wont be prime
            elif abs(word[1]) == remaining_slices -1 and word[1] > 0: 
                for i in range(ends-2):
                    if ValidCup(word, i):
                        newWord = ["cup" + str(i) + "." + word[0], word[1]-1]
                        nextListofWords.append(newWord)
            else:
                for i in range(ends):
                    newWord = ["cap" + str(i) + "." + word[0], word[1]+1]
                    nextListofWords.append(newWord)
                for i in range(ends-2):
                    if ValidPos(word, i):
                        newWord = ["pos" + str(i) + "." + word[0], word[1]]
                        nextListofWords.append(newWord)
                    if ValidNeg(word, i):
                        newWord = ["neg" + str(i) + "." + word[0], word[1]]
                        nextListofWords.append(newWord)
                    if ValidCup(word, i): 
                        newWord = ["cup" + str(i) + "." + word[0], word[1]-1]
                        nextListofWords.append(newWord)
        remaining_slices -= 1
        currentListofWords = nextListofWords[:]
        
    ListofWords = []
    for word in currentListofWords:
        ends = 5 + 2*word[1]
        for i in range(ends-2):
            if ValidCup(word, i): 
                newWord = "cup" + str(i) + "." + word[0]
                nextListofWords.append(newWord)

    return ListofWords
    