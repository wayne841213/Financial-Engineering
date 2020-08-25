# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:42:13 2019

@author: 郭威漢
"""
import numpy as np

S0=50
K1=50
r=0.1
q=0
sigma = 0.4
T=0.25
nr=100
exercise_type='A' # A:American E:European

t=0
S_max_0=50

T=T-t
dt=T/nr


#%% Monte Carlo

mu= np.log(S0)+(r-q-0.5*sigma**2)*dt

def payoff(ST,S_max):
    """
    if ST>=K1 and option_type=='C':
        po= ST-K1
    """
    if S_max>ST :
        po= S_max-ST
    else:
        po= 0
    return(po)


Value=[]


for x in range(20):
    lnST_list = np.random.normal(np.log(S0)+(r-q-0.5*sigma**2)*dt, sigma*dt**0.5, 10000) # mean and standard deviation
    PayOff=[]
    totalPayoff=0
    for lnST in lnST_list:
        S_max_t=S_max_0
        for nt in range(nr-1):
            lnST = np.random.normal(lnST+(r-q-0.5*sigma**2)*dt, sigma*dt**0.5, 1)[0]
            ST=np.exp(lnST)
            S_max_t=max(ST,S_max_t)
        PayOff.append(payoff(ST,S_max_t))
        totalPayoff+=payoff(ST,S_max_t)
    
    expectedPayoff=totalPayoff/len(lnST_list)
    value=np.exp(-r*T)*expectedPayoff
    Value.append(value)

mean=np.mean(Value)
standard_deviation=np.std(Value)

print('------------------------------------------Monte Carlo')
print(str(mean))

print('Standard Deviation = '+str(standard_deviation))

percentRange= [(mean-2*standard_deviation),(mean+2*standard_deviation)]
print('95% Range = '+str(percentRange))

#percentRange= [(mean-1*standard_deviation),(mean+1*standard_deviation)]
#print('68% Range = '+str(percentRange))

#%% CRR binomial tree


u= np.exp(sigma*dt**(0.5))
d= np.exp(-sigma*dt**(0.5))

#u= np.exp((r-0.5*sigma**2)*dt+sigma*dt**(0.5))
#d= np.exp((r-0.5*sigma**2)*dt-sigma*dt**(0.5))

#u= 0.5*sigma**2*dt+1+sigma*dt**(0.5)
#d= 0.5*sigma**2*dt+1-sigma*dt**(0.5)

p= (np.exp((r-q)*dt)-d)/(u-d)

def S(n,j):
    return(S0*u**(n-j)*d**(j))

# Basic
    
# Tree S(n,j)
tree=[[S0]]
CV=[[[]]]
S_max_matrix=[[[S_max_0]]]


for nn in range(1,nr+1):
    col=[]
    col2=[]
    S_max_column=[]
    for j in range(nn+1):
        St=S(nn,j)
        col.append(St)
        col2.append([])
        S_max_list=[]
        if j==0:
            S_max_parents = S_max_matrix[nn-1][j]
        elif j==nn:
            S_max_parents = S_max_matrix[nn-1][j-1]
        else:
            S_max_parents = S_max_matrix[nn-1][j-1]+S_max_matrix[nn-1][j]

        for S_max_t in S_max_parents:
            S_max_t = round(S_max_t,4)
            St = round(St,4)
            if S_max_t >= St and S_max_t not in S_max_list:
                S_max_list.append(S_max_t)
            elif S_max_t < St and St not in S_max_list:
                S_max_list.append(St)
        S_max_column.append(S_max_list)
        
    tree.append(col)
    CV.append(col2)
    S_max_matrix.append(S_max_column)
    

# payoff at maturity
CV[-1]=[]

for SN,S_max_list in zip(tree[-1],S_max_column):
    cv=[]
    for S_max in S_max_list:
        cv.append(payoff(SN,S_max))
    CV[-1].append(cv)


# backward reduction

for n in range(nr-1,-1,-1):
    for j in range(n+1):
        for S_max in S_max_matrix[n][j]:
            if S_max < min(S_max_matrix[n+1][j]):
                X=-1
            else:
                X = S_max_matrix[n+1][j].index(S_max)
            Y = S_max_matrix[n+1][j+1].index(S_max)
            cv = np.exp(-r*dt)*(p*CV[n+1][j][X]+(1-p)*CV[n+1][j+1][Y])
            
            if exercise_type=='A':
                cv=max(cv,payoff(tree[n][j],S_max))
            if cv not in CV[n][j]:
                CV[n][j].append(cv)


print('------------------------------CRR binomial tree Basic')
print(CV[0][0][0])

#  BONUS1
# Tree S(n,j)
tree=[[S0]]
CV=[[[]]]

def Smax(n,j):
    S_max_list=[]
    ST=S(n,j)
    for nn in range(n-j,-1,-1):
        S_max=S0*u**nn
        if S_max >= ST:
            S_max_list.append(S_max)
    return(S_max_list)

for nn in range(1,nr+1):
    col=[]
    col2=[]
    for j in range(nn+1):
        St=S(nn,j)
        col.append(St)
        col2.append([])
        
    tree.append(col)
    CV.append(col2)

# payoff at maturity
CV[-1]=[]

for SN,x in zip(tree[-1],range(len(tree[-1]))):
    cv=[]
    for S_max in Smax(nr,x):
        cv.append(payoff(SN,S_max))
    CV[-1].append(cv)

# backward reduction
for n in range(nr-1,-1,-1):
    for j in range(n+1):
        for S_max in Smax(n,j):
            if S_max < min(Smax(n+1,j)):
                X = -1
            else:
                X = Smax(n+1,j).index(S_max)

            Y  = Smax(n+1,j+1).index(S_max)
            
            cv = np.exp(-r*dt)*(p*CV[n+1][j][X]+(1-p)*CV[n+1][j+1][Y])
            
            if exercise_type=='A':
                cv =max(cv,payoff(tree[n][j],S_max))
            if cv not in CV[n][j]:
                CV[n][j].append(cv)

print('-----------------------------CRR binomial tree Bonus1')
print(CV[0][0][0])

#%%  BONUS2
# Tree S(n,j)
tree=[[0]]
CV=[[[]]]

u= np.exp(sigma*dt**(0.5))
d= np.exp(-sigma*dt**(0.5))

#u= np.exp((r-0.5*sigma**2)*dt+sigma*dt**(0.5))
#d= np.exp((r-0.5*sigma**2)*dt-sigma*dt**(0.5))

#u= 0.5*sigma**2*dt+1+sigma*dt**(0.5)
#d= 0.5*sigma**2*dt+1-sigma*dt**(0.5)



#p= (np.exp((r)*dt)*u-1) / (np.exp((r)*dt)*(u-d))
p= (np.exp((r-q)*dt)-d)/(u-d)

pu=p*u
pd=(1-p)*d
for nn in range(1,nr+1):
    col=[]
    col2=[]
    for j in range(nn+1):
        col.insert(0, max(u**j-1,0))
        col2.append([])
        
    tree.append(col)
    CV.append(col2)

# payoff at maturity
A=[]

for power in range(len(tree[-1])):
    A.insert(0, max(u**power-1,0))

CV[-1]=A

# backward reduction

for n in range(nr-1,-1,-1):
    for j in range(n+1):
        k=j
        if j == n:
            k=j-1
        CV[n][j]=np.exp(-r*dt)*(pd*CV[n+1][j]+pu*CV[n+1][k+2])
        
        if exercise_type=='A':
            CV[n][j] =max(CV[n][j],tree[n][j])
        
        
        


print('-----------------------------CRR binomial tree Bonus2')
print(CV[0][0]*S0)

