# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 14:17:47 2019

@author: 郭威漢
"""

import numpy as np

S0=50
K=50
r=0.1
q=0
sigma = 0.4
T=0.25

mr=200
nr=100
exercise_type='A' # A:American E:European

S_max=100
S_min=0

T=T
dt=T/nr

#%% Implicit 

f= [[0 for x in range(mr+1)] for y in range(nr+1)]

dS= (S_max-S_min)/mr

f[-1]=[max(dS*x-K,0) for x in range(mr+1)]


#i rwo j column

def a(j):
    a= (r-q)*j*dt/2-(sigma*j)**2*dt/2
    return a

def b(j):
    b= 1+(sigma*j)**2*dt+r*dt
    return b

def c(j):
    c= -(r-q)*j*dt/2-(sigma*j)**2*dt/2
    return c

A=[[0 for x in range(mr+1)] for y in range(mr-1)]

for i in range(mr-1):
    for j in range(mr+1):
        if i==j:
            A[i][j]=a(i+1)
            A[i][j+1]=b(i+1)
            A[i][j+2]=c(i+1)

for x in A:
    x.reverse()
A.reverse()
    
A=np.matrix(A)

A=np.delete(A,0,1)
A=np.delete(A,-1,1)

#%% Explicit 










