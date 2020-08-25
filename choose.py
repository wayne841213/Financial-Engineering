# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 23:35:11 2019

@author: Wayne
"""

from fractions import Fraction

def choose(n,k):
    if k > n//2: k = n - k
    p = Fraction(1)
    for i in range(1,k+1):
        p *= Fraction(n - i + 1, i)
    return int(p)


choose(100,50)