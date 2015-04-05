#HW6, problem 5
import numpy      
import matplotlib.pyplot as plt     
import copy 

#Reusing code from HW4+HW5
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
    
    #slight change, don't need k
def GaussSeidel(A,b,er):
    d=len(A)
    x=numpy.zeros(d)
    xtemp=copy.copy(x)
    abser=2*er
    while abser>er:
        xtemp=copy.copy(x) 
        for i in range(d):
            x[i]= (1.0/A[i][i])*(b[i]-sumLower(A,x,i)-sumUpper(A, xtemp,i))        
        abser=numpy.max(numpy.abs(xtemp-x))
    return x

#Slight change to this too, don't need b
def findiff(h,a,sigma_a,D):
    gpoints=2*a/h
    m=numpy.zeros((gpoints,gpoints))
    for i in range(len(m)):
        m[i][i]=2*D/(h**2)+sigma_a
        if i<len(m)-1:
            m[i][i+1]=-1*D/(h**2)
            m[i+1][i]=-1*D/(h**2)
            
    print('Matrix initialized')
    return m

#Power iteration
def powerit(phi0,k0,A, Q0, simga_f, e1, e2):
    m=1
    phi_m=numpy.float_(phi0)/numpy.linalg.norm(phi0)
    k_m=k0
    Q_m=Q0
    ek=2*e1
    ephi=2*e2
    while(ephi>e2 or ek>e1):
        phi_m1=copy.copy(phi_m)  
        Q_m1=copy.copy(Q_m)
        k_m1=k_m
        Q_m = sigma_f*phi_m1 #Doesn't matter if I use phi_m or phi_m1, at this point they're the same
        k_m = k_m1*(numpy.sum(Q_m)/numpy.sum(Q_m1))       
        phi_m=GaussSeidel(A,Q_m1/k_m1,e2)
        
        ek=numpy.abs((k_m-k_m1)/k_m)
        ephi=numpy.max(numpy.abs((phi_m-phi_m1)/phi_m))       
        
        m+=1
    phi_n=phi_m/numpy.linalg.norm(phi_m)
    return k_m,phi_n,m

a=4
D=1
sigma_a=0.7
sigma_f=0.6
h=0.1
    
A=findiff(h,a,sigma_a,D)
phi0=numpy.ones(len(A), dtype=float)
Q0=sigma_f*copy.copy(phi0)
k0=0.1 #How to choose this? Random guess
k, phi, m = powerit(phi0, k0, A, Q0, sigma_f, 1*10**-4, 1*10**-4)

print "Eigenvalue, power iteration: ", k
print "Number of power iterations: ", m

F=sigma_f*numpy.matrix(numpy.identity(len(A)))
B=numpy.linalg.inv(A)*F

kvals, kvec = numpy.linalg.eig(numpy.asarray(B))
k2=numpy.max(kvals)
print "Eigenvalue, built in solvers: ", k2
i=numpy.argmax(kvals)
vec=kvec[i]

x=numpy.linspace(-a,a, num=2*a/h)

plt.figure()
plt.plot(x,phi)
plt.plot(x,vec, 'r')
plt.xlabel(r'$x$ [cm]')
plt.ylabel(r'$\phi(x)$')
plt.title('HW6, problem 6: Eigenvalues of diffusion eq. by power iteration')

plt.legend(['Power iteration', 'Built in solvers'])
