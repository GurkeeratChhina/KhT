#!/usr/bin/env python3
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

# load libraries
import math
import timeit
from time import time
from fractions import Fraction # for rational numbers as field of coefficients, not fully supported yet

# useful general functions that have nothing to do with this project
def find_first_index(mylist,function):
    """compute the index of the first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if function(v))

def find_first(mylist,function):
    """compute first element in the list mylist on which function is true. 
    Throws error if there is no such element"""
    return next(v for v in mylist if function(v))

def indexQ(mylist,val):
    """compute the index of the first element in the list mylist which is equal to val. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if v==val)

def indexMemberQ(mylist,val):
    """compute the index of the first element in the list mylist which contains val. 
    Throws error if there is no such element"""
    return next(i for i,v in enumerate(mylist) if val in v)

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

