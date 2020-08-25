# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:27:42 2019

@author: Wayne
"""

import numpy as np

S0=50
K1=55
r=0.1
q=0.03
T=0.5
P=6.5

nr=100
convergence_criterion=10e-8

option_type='P' # C:call P:put
exercise_type='A' # A:American E:European
model='binomial' # 'Black_Scholes'  'binomia'


import scipy.stats as st
#st.norm.ppf(.95)
def N(dType,K,sigma):
    
    if dType==1:
        d=(np.log(S0/K)+(r-q+0.5*sigma**2)*T) / (sigma*T**(0.5))
    if dType==2:
        d=(np.log(S0/K)+(r-q-0.5*sigma**2)*T) / (sigma*T**(0.5))  
        
    return(st.norm.cdf(d))
    
	
def bisection(model):
    a=0.0001
    b=1.5
    c=(a+b)/2.0
    while (b-a)/2.0 > convergence_criterion and model=='binomial':
        if (binomial(c)-P) == 0:
            return c
        elif (binomial(a)-P)*(binomial(c)-P) < 0:
            b = c
        else :
            a = c
        c = (a+b)/2.0

    while (b-a)/2.0 > convergence_criterion and model=='Black_Scholes':
        if (Black_Scholes(c)-P) == 0:
            return c
        elif (Black_Scholes(a)-P)*(Black_Scholes(c)-P) < 0:
            b = c
        else :
            a = c
        c = (a+b)/2.0
    if model != 'binomial' and model != 'Black_Scholes':
        print('Error!')
    return(c)

def Black_Scholes(sigma):
    if option_type=='C':
        P= S0*np.exp(-q*T)*N(1,K1,sigma) -K1*np.exp(-r*T)*N(2,K1,sigma)
    elif option_type=='P':
        P= K1*np.exp(-r*T)*(1-N(2,K1,sigma)) -S0*np.exp(-q*T)*(1-N(1,K1,sigma)) 
    else:
        P = 0
        print('Error!')
    return(P)
    



def payoff(ST,option_type):
    if ST>=K1 and option_type=='C':
        po= ST-K1
    elif ST<=K1 and option_type=='P':
        po= K1-ST 
    else:
        po=0
    return(po)


# CRR binomial tree
dt=T/nr

def binomial(sigma):
    u= np.exp(sigma*dt**(0.5))
    d= np.exp(-sigma*dt**(0.5))
    p= (np.exp((r-q)*dt)-d)/(u-d)
    def S(n,j):
        return(S0*u**(n-j)*d**(j))
    # Tree S(n,j)
    tree=[]
    CV=[]
    for nn in range(nr+1):
        col=[]
        col2=[]
        for j in range(nn+1):
            col.append(S(nn,j))
            col2.append(0)
        tree.append(col)
        CV.append(col2)

    # payoff at maturity
    CV[-1]=[]
    
    for SN in tree[-1]:
        CV[-1].append(payoff(SN,option_type))
    
    # backward reduction
    for n in range(nr-1,-1,-1):
        for j in range(n+1):
            CV[n][j]=np.exp(-r*dt)*(p*CV[n+1][j]+(1-p)*CV[n+1][j+1])
            if exercise_type=='A':
                CV[n][j]=max(CV[n][j],payoff(tree[n][j],option_type))
    return(CV[0][0])
    
print('bisection-------------------------')
print('Black_Scholes   '+str(bisection('Black_Scholes')))
print('binomial tree   '+str(bisection('binomial')))

import math

def vega(sigma):
    
    d=(np.log(S0/K1)+(r-q+0.5*sigma**2)*T) / (sigma*T**(0.5))
    v= S0*(np.exp(-0.5*d**2)/((2*math.pi)**0.5))*T**0.5

    return(v)

def d_binomial(sigma):
    h= 10e-9
    d= ( binomial(sigma+h)-binomial(sigma))/h
    
    return(d)

def newton(model):
    x0= 1.5
    xn = x0
    max_iter=100000
    if model=='Black_Scholes':
        f=Black_Scholes
        Df=vega
    elif model=='binomial':
        f=binomial
        Df=d_binomial
    for n in range(0,max_iter):
        fxn = f(xn)-P
        if abs(fxn) < convergence_criterion:
            #print('Found solution after',n,'iterations.')
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            print('Zero derivative. No solution found.')  
            return None
        xn = xn - fxn/Dfxn
        if xn<0:
            xn=0.5
    print('Exceeded maximum iterations. No solution found.')
    return None
print('newton----------------------------')
print('Black_Scholes   '+str(newton('Black_Scholes')))
print('binomial tree   '+str(newton('binomial')))