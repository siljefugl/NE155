
#HW6, problem 4
import numpy      
import matplotlib.pyplot as plt     
import matplotlib.axes as ax     

import copy        

#Reusing code from HW5, but changed it a bit - the old code was quite ineffective
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
    
def Jacobi(A, b, er):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=copy.copy(x)
    k=0
    abser=2*er
    while abser>er:
        xtemp=copy.copy(x)
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,xtemp,i)-sumUpper(A, xtemp,i))
        
        abser=numpy.max(numpy.abs(xtemp-x))
        k+=1
    return x,k-1
    
def GaussSeidel(A,b,er):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=copy.copy(x)
    k=0
    abser=2*er
    while abser>er:
        xtemp=copy.copy(x) 
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))        
        abser=numpy.max(numpy.abs(xtemp-x))
        k+=1
    return x,k-1

def SOR(A,b,er,w):
    d=len(A)
    xtemp=numpy.zeros(d)
    x=numpy.zeros(d)
    k=0
    abser=2*er
    while abser>er:
        xtemp=copy.copy(x)
        for i in range(d):
            x[i]= (1-w)*xtemp[i]+ (w/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))             
        k+=1
        abser=numpy.max(numpy.abs(xtemp-x))

    return x, k-1

a=4
D=1
sigma=0.2
S=8
def fanalytical(a,sigma,S,x):
    return ((S/sigma)*(1-(numpy.exp(numpy.sqrt(sigma)*x)+numpy.exp(-numpy.sqrt(sigma)*x))/(numpy.exp(numpy.sqrt(sigma)*a)+numpy.exp(-numpy.sqrt(sigma)*a))))


#Note that this is not the general case of finite difference method, but a special case for our specific equation
def findiff(h,a,sigma,S):
    gpoints=2*a/h
    m=numpy.zeros((gpoints,gpoints))
    b=numpy.zeros(len(m))
    for i in range(len(m)):
        m[i][i]=2+(h**2)*sigma/D
        b[i]=(h**2)*S/D
        if i<len(m)-1:
            m[i][i+1]=-1
            m[i+1][i]=-1 
    print('Initalized vectors')
    return m,b
            
h=1
absolute_error=1*10**-5
print 'For h=',h
print '\n'
A,b=findiff(h,a,sigma,S)    
x=numpy.linspace(-a,a,num=len(b))
#sol=fanalytical(a,sigma,S,x) 

#4a) Jacobi method
SOLJacobi,kJacobi=Jacobi(A,b,absolute_error)
print 'Jacobi method; # of iterations:', kJacobi 
print  '\n' 

#4b) Gauss Seidel method
SOLGS,kGS=GaussSeidel(A,b,absolute_error)
print 'Gauss Seidel method; # of iterations:', kGS 
print  '\n' 

#4c) SOR method with omega=1.2
SOLSOR,kSOR=SOR(A,b,absolute_error, 1.2)
print 'SOR method; # of iterations:', kSOR
print  '\n'

hval=[1,0.5,0.1,0.05]
itJacobi_e1=[54, 171, 2243, 5777]
itGS_e1=[30, 94, 1320, 3676]
itSOR_e1=[19, 65, 957, 2758]
itJacobi_e2=[83, 281, 4872, 16235]
itGS_e2=[45, 149, 2635, 8905]
itSOR_e2=[28, 101, 1832, 6243]


plt.figure()
plt.title('HW 6, problem 4')
plt.xlabel('h-value [cm]')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('# of iterations needed')
plt.plot(hval,itJacobi_e1, 'ob')
plt.plot(hval,itGS_e1, 'or')
plt.plot(hval,itSOR_e1, 'og')
plt.plot(hval,itJacobi_e2, 'xb')
plt.plot(hval,itGS_e2, 'xr')
plt.plot(hval,itSOR_e2, 'xg')

plt.legend(['Jacobi e1','Gauss-Seidel e1','SOR e1','Jacobi e2','Gauss-Seidel e2','SOR e2']) 

