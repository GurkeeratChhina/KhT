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

from CobComplexes import *
from BNComplexes import *
from Drawing import *

def StringToStringList(str):
    list=[[word[0:3],int(word[3:])] for word in str.split('.')]
    list.reverse()
    return list
    
def StringListToString(strlst):
    strlst.reverse()
    newstring = ""
    for word in strlst[:-1]:
        newstring += word[0] + str(word[1]) + "."
    newstring += strlst[-1][0] + str(strlst[-1][1])
    return newstring

class Tangle(object): #TODO: implement orientations for all methods
    """ A tangle is either a string that goes through all the slices, seperated by '.' or a list of slices
        'slices' is a string with the left most slice being the bottom of the tangle
        'stringlist' is a list with the left most slice being the top of the tangle
        both formats include an integer as the index at which the slice is being added
        'pos'/'neg' count the number of positively/negatively oriented crossings, respectively
        note that pos and neg are also used to specify the unoriented crossing in the slices
        such as pos1 meaning the overstrand has positive slope, at index 1
        'pos' and 'neg' are calculated automatically from the given input orientation
        the input_orientations is a length n+1 list, where n is the number of tangle ends
        the first n entries specify the orientations on the n tangle ends, read left to right, top to bottom
        the last entry is a list of orientations on any closed components, in the order they are encountered top to bottom
        if there are no closed components, then the last element is an empty list
        an orientation is an element of [-1, 0, 1], where -1 means downward, +1 means upwards, and 0 means unoriented
        on closed components, a -1 means left at the first cap encountered, and a +1 means right
        'orientations' stores the orientation information at each of the slices, as a list of lists of orientations
        'orientations' is of length 'height'+1, where 'height' is the number of slices. 
        
        A tangle is oriented by calling the OrientTangle method
    """
    __slots__ = 'slices', 'stringlist', 'pos', 'neg', 'top', 'bot', 'orientations', 'height'
    
    
    def __init__(self, string = None, strlist = None, topends = 1, botends =3, orientation = None, pos = 0, neg = 0):
        if string is not None and strlist is not None: #initializing using both string and stringlist
            self.slices = string
            self.stringlist = strlist
        if string is not None: #initializing using only string
            self.slices = string
            self.stringlist = StringToStringList(string)
        elif strlist is not None: #initializing using only stringlist
            self.stringlist = strlist
            self.slices = StringListToString(strlist)
        else: # no data to initialize with
            raise Exception('Constructing an empty tangle')
        
        self.top = topends
        self.bot = botends
        self.height = len(self.stringlist)
        
        if orientation is None: # initializing orienation data with 0's for unoriented
            ends = topends
            SliceOrientations = [[0]*ends]
            for slice in self.stringlist:
                if slice[0] == "cup":
                    ends -= 2
                elif slice[0] == "cap":
                    ends += 2
                SliceOrientations += [[0]*ends]
            self.orientations = SliceOrientations
            self.pos = 0
            self.neg = 0
        else:
            self.orientations = orientation
            self.pos = pos
            self.neg = neg
    
    @classmethod
    def PretzelTangle(cls, left_twists, right_twists):
        tangle = "cup1."
        for i in range(abs(left_twists)):
            if left_twists < 0:
                tangle += "neg0."
            if left_twists > 0:
                tangle += "pos0."
        for i in range(abs(right_twists)):
            if right_twists < 0:
                tangle += "neg2."
            if right_twists > 0:
                tangle += "pos2."
        tangle += "cap3.cap1"
        return Tangle(tangle)
        
    @classmethod
    def LiamsTangle(cls, N, Twists):
        tangle = "cup1.neg0.pos2."
        for j in range(N):
            tangle += "neg0.neg1."
            for k in range(abs(Twists[j])):
                if Twists[j] < 0:
                    tangle += "neg2."
                if Twists[j] >= 0:
                    tangle += "pos2."
            tangle += "neg1.neg0."
        tangle += "pos2.cap3.neg0.neg0.cap1"
        return Tangle(tangle)
        
    @classmethod
    def quotient_of_2_m3_pretzel_tangle(cls, N):
        tangle = "cup2.cup3.neg1.neg2.neg0.neg1.cup4.pos3.pos2.neg4.neg3.cap5.cap4.pos1.cup3.cup4.pos2.pos3.neg4.neg5.cap6."
        for j in range(abs(N)):
            if N > 0:
                tangle += "pos3.pos4.pos0.pos3."
            if N < 0:
                tangle += "neg3.neg4.neg0.neg3."
        tangle += "cap1.cap0.cap1"
        return Tangle(tangle)
    
    @classmethod
    def two_twist_hitch(cls, m):
        tangle = "cup1.cup1.pos4.neg0.neg3.pos1.cap2."
        for j in range(abs(m)):
            if m > 0:
                tangle += "pos1."
            if m < 0:
                tangle += "neg1."
        tangle += "pos2.pos2.neg0.neg0.cap3.cap1"
        return Tangle(tangle)

    # ONLY USE THIS ON UNORIENTED TANGLES
    def OrientTangle(self, input_orientations):    
        if len(input_orientations) != self.top + self.bot +1:
            raise Exception("Number of specified orientations does not agree with number of tangle ends.")
        for index, boundary_orientation in enumerate(input_orientations[:-1]):
            if index < self.top:  #Propogate orientations from the top
                if self.orientations[0][index] == boundary_orientation: #strand already oriented
                    continue
                elif self.orientations[0][index] == 0: # unoriented strand to start
                    self.PropogateOrientations(1, boundary_orientation, 0, index)
                else: #opposite orientation there
                    raise Exception("Orientations on boundary are not consistent")
            else: #Propogate orientations from the bot
                if self.orientations[-1][index-self.top] == boundary_orientation: #strand already oriented
                    continue
                elif self.orientations[-1][index-self.top] == 0: #strand unoriented to start
                    self.PropogateOrientations(-1, boundary_orientation, self.height, index-self.top)
                else: #opposite orientation there
                    raise Exception("Orientations on boundary are not consistent")
        for closed_orientation in input_orientations[-1]:
            for index, oriented_slice in enumerate(self.orientations):
                if contains_0(oriented_slice):
                    self.PropogateOrientations(1, closed_orientation, index, indexQ(oriented_slice, 0))
                    break
        for x, slice in enumerate(self.stringlist):
            if slice[0] == "pos":
                left_end = self.orientations[x][slice[1]]
                right_end = self.orientations[x][slice[1]+1]
                if left_end*right_end == 1:
                    self.pos += 1
                else:
                    self.neg += 1
            elif slice[0] == "neg":
                left_end = self.orientations[x][slice[1]]
                right_end = self.orientations[x][slice[1]+1]
                if left_end*right_end == 1:
                    self.neg += 1
                else:
                    self.pos += 1
        
    def PropogateOrientations(self, initial_direction, initial_orientation, initial_height, initial_index):
        direction = initial_direction # either +1 or -1, meaning travelling down (resp. up) the tangle
        curr_orientation = initial_orientation
        max_height = self.height
        curr_height = initial_height
        curr_index = initial_index
        while True:
            self.orientations[curr_height][curr_index] = curr_orientation
            if direction == 1: # Move down strand, changing orientation, height, index, and direction if needed
                next_slice = self.stringlist[curr_height] #find the next slice to calculate new stuff
                if next_slice[0] == "cap" and next_slice[1] <= curr_index: 
                    curr_height += 1
                    curr_index += 2
                elif next_slice[0] == "cap" and next_slice[1] > curr_index:
                    curr_height += 1
                elif next_slice[0] == "pos" and next_slice[1] == curr_index:
                    curr_height += 1
                    curr_index += 1
                elif next_slice[0] == "pos" and next_slice[1] == curr_index-1:
                    curr_height += 1
                    curr_index -= 1
                elif next_slice[0] == "neg" and next_slice[1] == curr_index:
                    curr_height += 1
                    curr_index += 1
                elif next_slice[0] == "neg" and next_slice[1] == curr_index-1:
                    curr_height += 1
                    curr_index -= 1
                elif next_slice[0] == "cup" and next_slice[1] > curr_index:
                    curr_height += 1
                elif next_slice[0] == "cup" and next_slice[1] < curr_index -1:
                    curr_height += 1
                    curr_index -= 2
                elif next_slice[0] == "cup" and next_slice[1] == curr_index:
                    curr_index +=1
                    direction *= -1
                    curr_orientation *= -1
                elif next_slice[0] == "cup" and next_slice[1] == curr_index -1:
                    curr_index -= 1
                    direction *= -1
                    curr_orientation *= -1
                else:
                    curr_height +=1
            elif direction == -1: # Move up strand, changing orientation, height, index, and direction if needed
                next_slice = self.stringlist[curr_height-1]
                if next_slice[0] == "cup" and next_slice[1] <= curr_index:
                    curr_height -= 1
                    curr_index += 2
                elif next_slice[0] == "cup" and next_slice[1] > curr_index:
                    curr_height -= 1
                elif next_slice[0] == "pos" and next_slice[1] == curr_index:
                    curr_height -= 1
                    curr_index += 1
                elif next_slice[0] == "pos" and next_slice[1] == curr_index-1:
                    curr_height -= 1
                    curr_index -= 1
                elif next_slice[0] == "neg" and next_slice[1] == curr_index:
                    curr_height -= 1
                    curr_index += 1
                elif next_slice[0] == "neg" and next_slice[1] == curr_index-1:
                    curr_height -= 1
                    curr_index -= 1
                ###### TODO : CHECK THE FOLLOWING
                elif next_slice[0] == "cap" and next_slice[1] > curr_index:
                    curr_height -= 1
                elif next_slice[0] == "cap" and next_slice[1] < curr_index-1:
                    curr_height -= 1
                    curr_index -=2
                elif next_slice[0] == "cap" and next_slice[1] == curr_index:
                    curr_index +=1
                    direction *= -1
                    curr_orientation *= -1
                elif next_slice[0] == "cap" and next_slice[1] == curr_index -1:
                    curr_index -= 1
                    direction *= -1
                    curr_orientation *= -1
                else:
                    curr_height -=1               
            if curr_height == 0 or curr_height == max_height:
                break
            if self.orientations[curr_height][curr_index] != 0:
                break 

    def shift(self, n): #shifts all the strands of a tangle over by n, to be used in the vertical and horizontal sums
        def shift_letter(letter):
            return [letter[0], letter[1] + n]
        shiftlist = [shift_letter(letter) for letter in self.stringlist]
        return Tangle(None, shiftlist, self.top, self.bot, self.orientations, self.pos, self.neg)

    def vertical_sum(self, other): # Assumes self and other are 1-3 tangles. Converts them to 2-2 tangles and stacks them vertically #TODO: orientations
        #TODO: verify that self and other are 1-3
        return Tangle("cup2." + other.slices + "." + self.slices)    
 
    def horizontal_sum(self, other): # Assumes self and other are 1-3 tangles. Converts them to 2-2 tangles and stacks them horizontally #TODO: orientations
        #TODO: verify that self and other are 1-3
        return Tangle("cup1." + other.shift(2).slices + "." + self.slices)
    
    def toReduced_BNComplex(self, max_iter = 100, field = 2, options = "unsafe", intermediate_cleanup = False):
        """ Computes the reduced BN complex from the tangle self
            max_iter is the maximum number of iterations in the cleanup procedure
            start specifies the number of ends at the top
            field is as usual: 0 = Q, 1 = Z, p = F_p
            options is either "safe" or "unsafe" - if "safe" then it checks if the complex is valid at each step
            intermediate_cleanup will convert to a BN complex, and apply the cleanup lemma if there is an
            intermediate point where the tangle is a 4-ended tangle
        """
        cx= BNbracket(self.slices, self.pos, self.neg, self.top)
        BN_complex= cx.ToBNAlgebra(field)
        BN_complex.clean_up(max_iter)
        return BN_complex
    
    def Cable(self): # TODO: orientations
        if self.top != 1 or self.bot != 1:
            raise Exception("Tangle is not a 1-1 tangle to cable")
        NewStringList = [["cap", 1]]
        for letter in self.stringlist:
            if letter[0] in ["pos", "neg"]:
                NewStringList.extend([ [letter[0], 2*letter[1]+1], [letter[0], 2*letter[1]], [letter[0], 2*letter[1]+2], [letter[0], 2*letter[1]+1] ])
            elif letter[0] == "cap":
                NewStringList.extend([ ["cap", 2*letter[1]], ["cap", 2*letter[1] + 1] ])
            elif letter[0] == "cup":
                NewStringList.extend([ ["cup", 2*letter[1]+1], ["cup", 2*letter[1]] ])    
        
        return Tangle(None, NewStringList)
        
    def draw(self, filename, style="plain"):
        drawtangle(self.slices,filename,style,self.top)
    
    
 
 
##########################################################################################################
###                                                                                                    ###
###    The following code is incomplete, it is used to generate all prime tangles of a given length    ###
###    Note this does not mean the number of crossings, but the number of horizontal slices            ###
###    The output should be as a list of tangle words or as a list of tangles                          ###
###                                                                                                    ###
##########################################################################################################

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
        
    ListofWords = [] #TODO: CHECK FOR NULL
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


