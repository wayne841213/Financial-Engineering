# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:27:42 2019

@author: Wayne
"""

import numpy as np

S0=50
K1=50
r=0.1
q=0.05
sigma = 0.4
T=0.5
nr=100
option_type='P' # C:call P:put
exercise_type='E' # A:American E:European



import scipy.stats as st
#st.norm.ppf(.95)
    
def N(dType,K):
    
    if dType==1:
        d=(np.log(S0/K)+(r-q+0.5*sigma**2)*T) / (sigma*T**(0.5))
    if dType==2:
        d=(np.log(S0/K)+(r-q-0.5*sigma**2)*T) / (sigma*T**(0.5))  
        
    return(st.norm.cdf(d))



#%% Black-Scholes
    
if option_type=='C':
    ana= S0*np.exp(-q*T)*N(1,K1) -K1*np.exp(-r*T)*N(2,K1)
    
elif option_type=='P':
    ana= K1*np.exp(-r*T)*(1-N(2,K1)) -S0*np.exp(-q*T)*(1-N(1,K1))

else:
    ana='error'
print('----------------------------------------Black-Scholes')
print(str(ana))
#%% Monte Carlo

mu= np.log(S0)+(r-q-0.5*sigma**2)*T

def payoff(ST,option_type):
    if ST>=K1 and option_type=='C':
        po= ST-K1
    elif ST<=K1 and option_type=='P':
        po= K1-ST 
    else:
        po=0
    return(po)


Value=[]


for x in range(20):
    lnST = np.random.normal(mu, sigma*T**0.5, 10000) # mean and standard deviation
    ST=np.exp(lnST)
    
    PayOff=[]
    totalPayoff=0
    for s in ST:
        PayOff.append(payoff(s,option_type))
        totalPayoff+=payoff(s,option_type)
    
    expectedPayoff=totalPayoff/len(ST)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)

mean=np.mean(Value)
standard_deviation=np.std(Value)

#print('------------------------------------------Monte Carlo')
#print(str(mean))

#print('Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
#print('95% Range = '+str(percentRange))

percentRange= [(mean-1*standard_deviation),(mean+1*standard_deviation)]
#print('68% Range = '+str(percentRange))

#%% CRR binomial tree

dt=T/nr

u= np.exp(sigma*dt**(0.5))
d= np.exp(-sigma*dt**(0.5))

#u= np.exp((r-0.5*sigma**2)*dt+sigma*dt**(0.5))
#d= np.exp((r-0.5*sigma**2)*dt-sigma*dt**(0.5))

#u= 0.5*sigma**2*dt+1+sigma*dt**(0.5)
#d= 0.5*sigma**2*dt+1-sigma*dt**(0.5)

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
            #print(str(CV[n][j])+' '+str(payoff(tree[n][j],option_type)))


print('------------------------------------CRR binomial tree')
print(CV[0][0])



#%% CRR one column Bonus1

CV1=[]

for SN in tree[-1]:
    CV1.append(payoff(SN,option_type))
    
for n in range(nr-1,-1,-1):
    for j in range(n+1):
        CV1[j]=np.exp(-r*dt)*(p*CV1[j]+(1-p)*CV1[j+1])
        if exercise_type=='A':
            CV1[j]=max(CV1[j],payoff(tree[n][j],option_type))
    CV1.pop(-1)
    
print(CV1[0])

#%% combinatorial method

from fractions import Fraction

def choose(n,k):
    if k > n//2: k = n - k
    p = Fraction(1)
    for i in range(1,k+1):
        p *= Fraction(n - i + 1, i)
    return p

nr=10000

cmV=0

aaa = Fraction(1)

for j in range(nr+1):
    if j != 0:
        aaa *= Fraction(nr - j +1 , j)*(1-p)/p
    if aaa*(p**(nr))>0 and aaa*(p**(nr))<1000:
        aaa *= (p**(nr))
        cmV *= (p**(nr))
    cmV+= aaa*max(S0*(u**(nr-j))*(d**j)-K1,0)
     
cmV=cmV*np.exp(-r*T)

print('---------------------------------combinatorial method')
print(cmV)



