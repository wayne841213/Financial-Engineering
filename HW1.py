# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:27:42 2019

@author: Wayne
"""

import numpy as np

S0=100
r=0.05
q=0.02
sigma = 0.5
T=0.4
K1=np.float64(90)
K2=np.float64(98)
K3=np.float64(102)
K4=np.float64(104)




import scipy.stats as st
#st.norm.ppf(.95)
    
def N(dType,K):
    
    if dType==1:
        d=(np.log(S0/K)+(r-q+0.5*sigma**2)*T) / (sigma*T**(0.5))
    if dType==2:
        d=(np.log(S0/K)+(r-q-0.5*sigma**2)*T) / (sigma*T**(0.5))  
        
    return(st.norm.cdf(d))



# Basic 2
    
ana = S0*np.exp(-q*T)*(N(1,K1)-N(1,K2))-K1*np.exp(-r*T)*(N(2,K1)-N(2,K2))+(K2-K1)*np.exp(-r*T)*(N(2,K2)-N(2,K3))+(K2-K1)/(K4-K3)*(K4*np.exp(-r*T)*(N(2,K3)-N(2,K4))-S0*np.exp(-q*T)*(N(1,K3)-N(1,K4)))

print('Analytical Solution = '+str(ana))


# Bonus

mu= np.log(S0)+(r-q-0.5*sigma**2)*T

def payoff(ST):
    if K2>=ST and ST>K1:
        po= ST-K1
    elif K3>=ST and ST>K2:
        po= K2-K1
    elif K4>=ST and ST>K3:
        po= (K4-ST)*(K2-K1)/(K4-K3)
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
        PayOff.append(payoff(s))
        totalPayoff+=payoff(s)
    
    expectedPayoff=totalPayoff/len(ST)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)


mean=np.mean(Value)
standard_deviation=np.std(Value)
print('Mean = '+str(mean)+'\n'+'Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
print('95% Range = '+str(percentRange))

percentRange= [(mean-1*standard_deviation),(mean+1*standard_deviation)]
print('68% Range = '+str(percentRange))


#%%Plot

#abs(mu - np.mean(s)) < 0.01
#abs(sigma - np.std(s, ddof=1)) < 0.01

sigma=sigma*T**0.5

import matplotlib.pyplot as plt 
count, bins, ignored = plt.hist(lnST, 3000, density=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
             np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
         linewidth=2, color='r')
plt.show()

mean=np.mean(lnST)
standard_deviation=np.std(lnST)
print('Mean = '+str(mean)+'\n'+'Standard Deviation = '+str(standard_deviation))
print('dleta Mean = '+str(abs(mean-mu))+'\n'+'delta Standard Deviation = '+str(abs(standard_deviation-sigma)))