# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 19:25:51 2015

@author: Silje
"""

import numpy     
from scipy import linalg              
import matplotlib.pyplot as plt

#Problem 4a)
A=numpy.zeros((100,100)) 
b=numpy.zeros((100,1))
A[0,0]=2
A[0,1]=-1

A[99,98]=-1
A[99,99]=2
b[99]=99

for j in range(1,99):
    A[j,j-1]=-1
    A[j,j]=2
    A[j,j+1]=-1
    b[j]=j
    
#4b)
c = numpy.linalg.cond(A)
print c

#4c)
Ainv=linalg.inv(A)
x1=Ainv.dot(b)

#4d)
x2=linalg.solve(A,b)

x=numpy.linspace(0,99,num=100)
plt.figure()
plt.plot(x1,'rx')
plt.plot(x2,'bx')
plt.title('Solutions of tridiagonal system, problem 4')
plt.xlabel('xth number of solution')
plt.ylabel('x-values')
plt.legend(['Inversion','linalg.solve'])
