#HW5, Problem 4

import numpy     
from scipy import linalg  
import copy        

def sumUpper(A,x, i):
    s=0
    if (i+1)<=len(A):
        for j in range(i+1, len(A)):
            s+=A[i][j]*x[j]      
    return s
    
def sumLower(A,x, i):
    s=0
    if ((i-1)>=0):
        for j in range(0, i):
            s+=A[i][j]*x[j] 
    return s
    
def Jacobi(A, b, k):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=copy.copy(x)
    l=0
    while l<k:
        xtemp=copy.copy(x)
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,xtemp,i)-sumUpper(A, xtemp,i))
        l+=1
    return x
    
def GaussSeidel(A,b,k):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=copy.copy(x)
    l=0
    while l<k:
        xtemp=copy.copy(x) 
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))        
        l+=1
    return x

def SOR(A,b,k,w):
    d=len(A)
    xtemp=numpy.zeros(d)
    x=numpy.zeros(d)
    l=0
    while l<k:
        xtemp=copy.copy(x)
        for i in range(d):
            x[i]= (1-w)*xtemp[i]+ (w/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))             
        l+=1
        
    return x

#Building matrix
n=5
eps=1*10**(-6)
A=numpy.zeros((n,n))
b=numpy.zeros(n)
for i in range(n-1):
    A[i][i]=4
    A[i][i+1]=-1
    A[i+1][i]=-1
    b[i]=100
A[n-1][n-1]=4
b[n-1]=100
sol=linalg.solve(A,b) 

#Jacobi method
er=10
rer=10
k=5
while rer>1*10**-6:
    solJ1=Jacobi(A,b,k)
    k+=1
    solJ2=Jacobi(A,b,k)
    er=numpy.amax(numpy.abs(sol-solJ2))
    rer=numpy.linalg.norm(solJ2-solJ1)/numpy.linalg.norm(solJ2)
    if k>100:
        er=0  
print 'With the Jacobi iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The relative error is:', rer, ', with k=', k
print 'The solution vector is:', solJ2
print  '\n' 

#Gauss Seidel method
er=10
rer=10
k=5
while rer>1*10**-6:
    solGS1=GaussSeidel(A,b,k)
    k+=1
    solGS2=GaussSeidel(A,b,k)
    er=numpy.amax(numpy.abs(sol-solGS2))
    rer=numpy.linalg.norm(solGS2-solGS1)/numpy.linalg.norm(solGS2)
    
    if k>100:
        er=0  
print 'With the Gauss Seidel iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The relative error is:', rer, ', with k=', k
print 'The solution vector is:', solGS2
print  '\n' 

#SOR method
rer=10
er=10
k=5
w=1.1
while rer>1*10**-6:
    solSOR1=SOR(A,b,k,w)
    k+=1
    solSOR2=SOR(A,b,k,w)
    er=numpy.amax(numpy.abs(sol-solSOR2))
    rer=numpy.linalg.norm(solSOR2-solSOR1)/numpy.linalg.norm(solSOR2)
    if k>100:
        rer=0  
print 'With the SOR iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The relative error is:', rer, ', with k=', k
print 'The solution vector is:', solSOR2

#SOR method, find best w (so error less than 10^-6 and k as small as possible)

omega=numpy.linspace(0.5,1.5,num=50)
kmin=8
wopt=1.1
eopt=7.6556794113*10**-06
for i in range(len(omega)):
    k=1
    while er>1*10**-6:
        solSOR1=SOR(A,b,k,omega[i])
        k+=1
        solSOR2=SOR(A,b,k,omega[i])
        er=numpy.amax(numpy.abs(sol-solSOR2))
        rer=numpy.linalg.norm(solSOR2-solSOR1)/numpy.linalg.norm(solSOR2)
        if k>100:
            rer=0 
    print k        
    if (k<kmin):
        print k
        kmin=copy.copy(k)
        eopt=copy.copy(er)
        wopt=omega[i]
    elif (k==kmin) and er<eopt:
        eopt=copy.copy(er)
        wopt=omega[i]


print 'The optimal omega is:', wopt

rer=10
er=10
k=1
while rer>1*10**-6:
    solSOR1=SOR(A,b,k,wopt)
    k+=1
    solSOR2=SOR(A,b,k,wopt)
    er=numpy.amax(numpy.abs(sol-solSOR2))
    rer=numpy.linalg.norm(solSOR2-solSOR1)/numpy.linalg.norm(solSOR2)
    if k>100:
        rer=0  

print 'The absolute error is now:', er, ', with k=', k
print 'The relative error is now:', rer, ', with k=', k
print 'The solution vector is now:', solSOR2

