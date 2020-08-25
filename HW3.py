# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:27:42 2019

@author: Wayne
"""

import numpy as np


K=100
r=0.1
T=0.5


nOfS=5

S0=[95]*nOfS
Q=[0.05]*nOfS
Sigma = [0.5]*nOfS





"""
C=[[1,-1],
   [-1,1]]
"""

C=[[1,0.5,0.5,0.5,0.5],
   [0.5,1,0.5,0.5,0.5],
   [0.5,0.5,1,0.5,0.5],
   [0.5,0.5,0.5,1,0.5],
   [0.5,0.5,0.5,0.5,1]
   ]


#% Cholesky Decomposition 
from math import sqrt
 
def cholesky(C):
    
    n = len(C)

    # Create zero matrix for L
    A = [[0.0] * n for i in range(n)]

    # Perform the Cholesky decomposition
    A[0][0]=sqrt(C[0][0])
    
    for j in range(1,n):
        A[0][j]=C[0][j]/A[0][0]
        
    for i in range(1,n):
        tmp_sum = sum(A[k][i]**2 for k in range(i))      
        A[i][i] = sqrt(C[i][i] - tmp_sum)
           
        for j in range(i+1,n):
            tmp_sum = sum(A[k][i] * A[k][j] for k in range(i)) 
            if i != j:
                A[i][j] = (1.0/A[i][i] * (C[i][j]-tmp_sum))
                
    A[n-1][n-1] = sqrt(C[n-1][n-1] - tmp_sum)
    return np.matrix(A)

A = cholesky(C)

#print( A.getT()*A )
#print( np.matrix(C) )


import scipy
import scipy.linalg   # SciPy Linear Algebra Library

#U = scipy.linalg.cholesky(C)
AA=np.matrix(U)


#% Monte Carlo



def payoff(ST):
    po=max(max(ST)-K,0)
    return(po)



Value=[]

#%%Basic
LNST=[]
for x in range(20):
    mu=[]
    lnST=[]
    for s0,q,sigma in zip(S0,Q,Sigma):
        mu.append( np.log(s0)+(r-q-0.5*sigma**2)*T )
        lnST.append( np.random.normal(0, sigma*T**0.5, 10000)) # mean and standard deviation
    lnST0= np.matrix(lnST).T
    
    LNST.append(lnST0)
    #print(lnST0.mean(0)) 

    lnST = lnST0 *A
    
    for n,lnst in zip(range(len(lnST)),lnST):
        lnst += np.matrix(mu)
        
    ST=np.exp(lnST)   
    
    PayOff=[]
    totalPayoff=0
    for s in ST:
        PayOff.append(payoff(s.A1))
        totalPayoff+=payoff(s.A1)
    expectedPayoff=totalPayoff/len(ST)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)


mean=np.mean(Value)
standard_deviation=np.std(Value)



print('------------------------------------------Basic')
print(str(mean))

print('Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
print('95% Range = '+str(percentRange))

percentRange= [(mean-1*standard_deviation),(mean+1*standard_deviation)]
#print('68% Range = '+str(percentRange))


#%%Bonus1

LNST1=[]

for x,lnST in zip(range(20),LNST):
    
    
    for n in range(int(len(lnST)/2)):
        lnST[len(lnST)-1-n] = -lnST[n]

    lnST1 =  lnST   / ( lnST.std(0) / (np.matrix(Sigma)*T**0.5) )
    LNST1.append(lnST1)
    
    lnST = lnST1* A
    
    for n,lnst in zip(range(len(lnST)),lnST):
        lnst += np.matrix(mu)
        
    ST=np.exp(lnST)   
    
    PayOff=[]
    totalPayoff=0
    for s in ST:
        PayOff.append(payoff(s.A1))
        totalPayoff+=payoff(s.A1)
    expectedPayoff=totalPayoff/len(ST)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)


mean=np.mean(Value)
standard_deviation=np.std(Value)



print('------------------------------------------Bonus1')
print(str(mean))

print('Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
print('95% Range = '+str(percentRange))


#%%Bonus2

for x,lnST1 in zip(range(20),LNST1):
    
    Cov = np.cov(lnST1.T)
    #U  = scipy.linalg.cholesky(Cov, lower=False)
    
    U =cholesky(Cov)
    AC =np.matrix(U)
    
    
    Var= np.matrix(Sigma).A1*T**0.5
    IVar = [[0.0] * len(Var) for i in range(len(Var))]
    for var,x in zip(Var,range(len(Var))):
        IVar[x][x]=var
    
    IVar=np.matrix(IVar)
    
    lnST = lnST1* AC.I * IVar *A
  
    
    for n,lnst in zip(range(len(lnST)),lnST):
        lnst += np.matrix(mu)
        
    ST=np.exp(lnST)   
    
    PayOff=[]
    totalPayoff=0
    for s in ST:
        PayOff.append(payoff(s.A1))
        totalPayoff+=payoff(s.A1)
    expectedPayoff=totalPayoff/len(ST)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)


mean=np.mean(Value)
standard_deviation=np.std(Value)



print('------------------------------------------Bonus2')
print(str(mean))

print('Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
print('95% Range = '+str(percentRange))