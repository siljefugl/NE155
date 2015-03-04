# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 20:31:09 2015

@author: Silje
"""
#HW4, Problem 6: Write a program to implement the Jacobi, Gauss Seidel, and 
#SOR methods for a matrix with n unknwons. Solve this system of equations with 
#each method using the program; use w=1.1 for SOR, x0=0 and n=5. Indicate the 
#number of iterations required to meet this tolerance for each method.

import numpy     
from scipy import linalg          

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
    xtemp=x
    l=0
    while l<k:
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,xtemp,i)-sumUpper(A, xtemp,i))
        xtemp=x 
        l+=1
    return x
    
def GaussSeidel(A,b,k):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=x
    l=0
    while l<k:
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))      
        xtemp=x   
        l+=1
        
    return x

def SOR(A,b,k,w):
    d=len(A)
    xtemp=numpy.zeros(d)
    x=numpy.zeros(d)
    l=0
    while l<k:
        for i in range(d):
            x[i]= (1-w)*xtemp[i]+ (w/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))          
        xtemp=x   
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
k=5
while er>1*10**-6:
    solJ=Jacobi(A,b,k)
    er=numpy.amax(numpy.abs(sol-solJ))
    k+=1
    if k>100:
        er=0  
print 'With the Jacobi iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The solution vector is:', solJ
print  '\n' 

#Gauss Seidel method
er=10
k=5
while er>1*10**-6:
    solGS=GaussSeidel(A,b,k)
    er=numpy.amax(numpy.abs(sol-solGS))
    k+=1
    if k>100:
        er=0  
print 'With the Gauss Seidel iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The solution vector is:', solGS
print  '\n' 

#SOR method
er=10
k=5
w=1.1
while er>1*10**-6:
    solSOR=SOR(A,b,k,w)
    er=numpy.amax(numpy.abs(sol-solSOR))
    k+=1
    if k>100:
        er=0  
print 'With the SOR iteration method:'
print 'The absolute error is:', er, ', with k=', k
print 'The solution vector is:', solSOR

