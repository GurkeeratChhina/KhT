#!/usr/bin/env python3

# load libraries
from itertools import groupby
from itertools import product
import itertools as itertools
import math

#TO DO: 
#4. Implement gaussian elimination
#5. Implement crossings
#6. Test all of the above features

# useful general functions that have nothing to do with this project
def find_first_index(mylist,function):
    """compute the index of the first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if function(v))

def find_first(mylist,function):
    """compute the index of the first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(v for v in mylist if function(v))

def indexQ(mylist,val):
    """compute the index of the first element in the list mylist which is equal to val. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if v==val)

def flatten(mylist):
    """analogue to Mathematica Flatten"""
    return [j for i in mylist for j in i]
    
def contains_0(list):
    """used to find basepoint containing component"""
    return 0 in list

def notcontains_0(list):
    """used to find the first component that doesn't contain the basepoint, usually for cobordisms with only 2 components"""
    return 0 not in list

#product, without using numpy
def prod(iterable):
    prod = 1
    for x in iterable:
        prod *= x
    return prod

#
# dictionary
#
# CLT= crossingless tangle
# TEI= tangle end index
# DS = alternate way of representing cobordisms of 4 ended tangles, using powers of D and powers of S
