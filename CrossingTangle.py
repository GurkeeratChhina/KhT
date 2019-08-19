import math
from KhT import *
from Tangles import *
from Complex import *
from Drawing import *
from Cobordisms import *
from BNComplexes import *

class Tangle(object):
    __slots__ = 'slices'
    
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
    
    def toReduced_BNComplex(self, max_iter = 100, start = 1, field = 2, options = "unsafe", intermediate_cleanup = False):
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
            
            if ends == 3 and intermediate_cleanup:
                BNcx = CobComplex2BNComplex(cx, field)
                BNcx.clean_up(max_iter)
                cx = BNComplex2CobComplex(BNcx)

                
            
        cx.shift_qhd(pos-2*neg,-neg,0.5*neg)
        print("Completed the computation successfully.                                              ")
        BN_complex= CobComplex2BNComplex(cx, field)
        BN_complex.clean_up(max_iter)
        return BN_complex
    
    def draw(self, filename, style="plain"):
        drawtangle(self.slices,filename,style)
    
    
    
def ValidCap(word, i): #check if adding cap at index i is valid
    ends = len(word[0])
    if i > ends - 2:
        return True
    previous_left = word[0][i]
    previous_right = word[0][i+1]
    if previous_left == ["cap", i] and previous_right == ["cap", i]:
        return False
    else:
        return True

def ValidCup(word, i): # check if adding cup at index i is valid
    previous_left = word[0][i]
    previous_right = word[0][i+1]
    if previous_left == ["cap", i-1]:
        return False
    elif previous_right == ["cap", i+1]:
        return False
    elif previous_left == ["cap", i] and previous_right == ["cap", i]:
        return False
    elif previous_left == ["neg", i] and previous_right == ["neg", i]:
        return False
    elif previous_left == ["pos", i] and previous_right == ["pos", i]:
        return False
    else:
        return True
        
def ValidPos(word, i): #check if adding pos at index i is valid
    previous_left = word[0][i]
    previous_right = word[0][i+1]
    if previous_left == ["neg", i] and previous_right == ["neg", i]:
        return False
    elif previous_left == ["cap", i] and previous_right == ["cap", i]:
        return False
    else:
        return True
        
def ValidNeg(word, i): #check if adding neg at index i is valid
    previous_left = word[0][i]
    previous_right = word[0][i+1]
    if previous_left == ["pos", i] and previous_right == ["pos", i]:
        return False
    elif previous_left == ["cap", i] and previous_right == ["cap", i]:
        return False
    else:
        return True


def GenerateTangleWords(max_length):
    """a Word is [ [strand info], [letter, index], [letter, index], ...  ]
        where [letter, index] describes the tangle word, read left to right
        and [strand info] is a list of [letter, index] such that [strand_info][i]
        gives the first thing the i-th strand encounters, from bottom up"""
    
    def shiftstrand(strand, amount):
        return [strand[0], strand[1] + amount]
    
    remaining_slices = max_length-1
    currentListofWords = [ [[["null", 0], ["cap", 1], ["cap", 1]], ["cap", 1]] ]
    while remaining_slices > 1:
        nextListofWords = []
        for word in currentListofWords:
            ends = len(word[0])
            if ends - 2*(remaining_slices - 1) > 5: # only adding cups wont get ends down enough
                continue #Do nothing
            elif ends == 1 and remaining_slices in [3,2]: # can only add caps to get enough ends in time
                continue #Do nothing, as wont be prime
            elif ends == 3 and remaining_slices == 2: # can only add 1 cap to get enough ends in time
                continue #Do nothing, as wont be prime
            elif ends - 2*(remaining_slices - 1) == 5: # can only add cups 
                for i in range(ends-2):
                    if ValidCup(word, i):
                        info = word[0] 
                        newinfo = info[:i] + [shiftstrand(strand, -2) for strand in info[i+2:]] #remove strand i, i+1, shifts everything after
                        newWord = [newinfo, ["cup", i]] + word[1:] #adds cupi to the word
                        nextListofWords.append(newWord)
            elif remaining_slices == 2: # only crossings
                for i in range(ends-2):
                    if ValidPos(word, i):
                        info = word[0] 
                        newinfo = info[:i] + [["pos", i], ["pos", i]] + info[i+2:] # changes strand i, i+1 to posi
                        newWord = [newinfo, ["pos", i]] + word[1:] #adds posi to the word
                        nextListofWords.append(newWord)
                    if ValidNeg(word, i):
                        info = word[0] 
                        newinfo = info[:i] + [["neg", i], ["neg", i]] + info[i+2:] # changes strand i, i+1 to negi
                        newWord = [newinfo, ["neg", i]] + word[1:] #adds negi to the word
                        nextListofWords.append(newWord)
            else: # Anything goes
                for i in range(ends):
                    if ValidCap(word, i):
                        info = word[0]
                        newinfo = info[:i] + [["cap", i], ["cap", i] ] + [shiftstrand(strand, 2) for strand in info[i:]] #adds cup at strand i, i+1, shifts everything after
                        newWord = [newinfo, ["cap", i]] + word[1:] #adds capi to the word
                        nextListofWords.append(newWord)
                for i in range(ends-2):
                    if ValidPos(word, i):
                        info = word[0] 
                        newinfo = info[:i] + [["pos", i], ["pos", i]] + info[i+2:] # changes strand i, i+1 to posi
                        newWord = [newinfo, ["pos", i]] + word[1:] #adds posi to the word
                        nextListofWords.append(newWord)
                    if ValidNeg(word, i):
                        info = word[0] 
                        newinfo = info[:i] + [["neg", i], ["neg", i]] + info[i+2:] # changes strand i, i+1 to negi
                        newWord = [newinfo, ["neg", i]] + word[1:] #adds negi to the word
                        nextListofWords.append(newWord)
                    if ValidCup(word, i): 
                        info = word[0]
                        newinfo = info[:i] + [shiftstrand(strand, -2) for strand in info[i+2:]] #remove strand i, i+1, shifts everything after
                        newWord = [newinfo, ["cup", i]] + word[1:] #adds cupi to the word
                        nextListofWords.append(newWord)

        remaining_slices -= 1
        currentListofWords = nextListofWords[:]
        
    ListofWords = [] #CHECK FOR NULL
    for word in currentListofWords:
        ends = len(word[0])
        for i in range(ends-2):
            if ValidCup(word, i): 
                info = word[0]
                newinfo = info[:i] + [shiftstrand(strand, -2) for strand in info[i+2:]] #remove strand i, i+1, shifts everything after
                newWord = [newinfo, ["cup", i]] + word[1:] #adds cupi to the word
                ListofWords.append(newWord)
    
    def primeword(word):
        info = word[0]
        for index, strand in enumerate(info):
            if strand[0] == "null":
                return False
            elif index == len(info) -1:
                continue
            elif strand[0] in ["cap", "pos", "neg"]:
                if strand == info[index+1]:
                    return False
        return True
        
    return [ word for word in ListofWords if primeword(word)]

def TangleWordToTangle(word):
    string = ""
    for letter in word[1:]:
        string += letter[0] + str(letter[1]) + "."
    string = string[:-1]
    return Tangle(string)