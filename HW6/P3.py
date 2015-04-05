#HW6, problem 3
import numpy      
import matplotlib.pyplot as plt
import math

a=4
D=1
A=0.2
S=8
def fanalytical(a,sigma,S,x):
    return ((S/sigma)*(1-(numpy.exp(numpy.sqrt(sigma)*x)+numpy.exp(-numpy.sqrt(sigma)*x))/(numpy.exp(numpy.sqrt(sigma)*a)+numpy.exp(-numpy.sqrt(sigma)*a))))

hval=[0.01,0.05,0.1,0.5,1]
relermax=numpy.zeros(len(hval))
relermin=numpy.zeros(len(hval))
nmeshes=numpy.zeros(len(hval))
for j in range(len(hval)):
    gpoints=2*a/hval[j]
    nmeshes[j]=gpoints
    m=numpy.zeros((gpoints,gpoints))
    b=numpy.zeros(len(m))
    for i in range(len(m)):
        m[i][i]=2+hval[j]**2*A/D
        b[i]=hval[j]**2*S/D
        if i<len(m)-1:
            m[i][i+1]=-1
            m[i+1][i]=-1       
    fnum=numpy.linalg.solve(m,b)
    x=numpy.linspace(-a,a,num=gpoints)
    fan=fanalytical(a,A,S,x) 
    reler=numpy.abs((fan-fnum)/fan)
    #But don't want outer points as they are infinity
    relermax[j]=numpy.max(reler[1:(len(reler)-1)])
    relermin[j]=numpy.min(reler)


plt.figure()
plt.plot(nmeshes,relermax, 'o')
plt.plot(nmeshes,relermin, 'or')
plt.title('HW 6, problem 3')
plt.xlabel('# of meshes')
plt.ylabel('Relative error')
plt.legend(['Max','Min']) 
