#HW6, problem 2
import numpy      
import matplotlib.pyplot as plt

a=4
D=1
A=0.2
S=8
h=0.1

gpoints=2*a/h
m=numpy.zeros((gpoints,gpoints))
b=numpy.zeros(len(m))
for i in range(len(m)):
    m[i][i]=2+h**2*A/D
    b[i]=h**2*S/D
    if i<len(m)-1:
        m[i][i+1]=-1
        m[i+1][i]=-1
        

f=numpy.linalg.solve(m,b)
x=numpy.linspace(-a,a,num=gpoints)

plt.figure()
plt.plot(x,f,'b')

def fanalytical(a,sigma,S,x):
    return ((S/sigma)*(1-(numpy.exp(numpy.sqrt(sigma)*x)+numpy.exp(-numpy.sqrt(sigma)*x))/(numpy.exp(numpy.sqrt(sigma)*a)+numpy.exp(-numpy.sqrt(sigma)*a))))

plt.plot(x,fanalytical(a,A,S,x),'r')
plt.title('HW 6, problem 2')
plt.xlabel('x')
plt.ylabel(r'$\phi$(x)')
plt.legend(['Numerical solution','Analytical solution'])    

