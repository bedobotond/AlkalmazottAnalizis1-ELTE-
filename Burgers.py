import scipy.optimize
from numpy import * 
from scipy import integrate
import numpy as np
import scipy
%matplotlib inline
import matplotlib.pyplot as plt
import mpld3
mpld3.enable_notebook()

# Problem: nonviscous 1D Burgers with periodic B.C.
# Goal: Investigation of CLAWPACK  and semianalytic solution for this problem

## Characteristic function
def func(xi):
    return (np.sin(2.0*np.pi*xi)+0.50)*(t-0)+xi-x

## Parameters
a=0.0
b=1.0
N=2**2
t=0.18
h=(b-a)/N

## Semianalytic solution with methods of characteristic
xx=(np.asarray(range(1,2*N,2)))/(2.0*N)
xi = zeros(len(xx))
sol = [] 
for guess in range(N):
    for i in range(guess):
        x = xx[i]
        xi[i] = scipy.optimize.fsolve(func, xi[i-1]) 
    
    xi[N-1]=xi[0]+1
    for j in np.arange(len(xx)-2, guess-1,-1):
        x = xx[j]
        xi[j] = scipy.optimize.fsolve(func, xi[j+1]) 
    q_0_xi=np.sin(2.0*np.pi*xi)+0.5
    sol.append(np.abs(0.50-scipy.integrate.trapz(q_0_xi,xx)))
(m,index1) = min((index2,index1) for index1,index2 in enumerate(sol))

##
h1=(b-a)/N
Semi_1_Before=q_0_xi

## CLAWPACK solution
from clawpack.pyclaw import examples
claw = examples.burgers_1d.setup(mx=N)     
claw.tfinal = 0.16#58214
claw.run()
Numeric_1_Before=claw.frames[10].q[0,:]

## Error definition
q=3.0
E_h=(h1*np.sum(np.abs(Numeric_1_Before-Semi_1_Before)**q))**(1.0/q)
print E_h

## Repeating the process on a finer grid
xx=(np.asarray(range(1,4*N,2)))/(4.0*N)
t=0.16#58214
xi = zeros(len(xx))
for i in range(guess):
    x = xx[i]
    xi[i] = scipy.optimize.fsolve(func, xi[i-1])
xi[2*N-1]=xi[0]+1
for j in np.arange(len(xx)-2, guess-1,-1):
    x = xx[j]
    xi[j] = scipy.optimize.fsolve(func, xi[j+1]) 
q_0_xi=np.sin(2.0*np.pi*xi)+0.5
h2=(b-a)/(2.0*N)
Semi_2_Before=q_0_xi
from clawpack.pyclaw import examples
claw = examples.burgers_1d.setup(mx=2*N)     
claw.tfinal = 0.16#58214
claw.run()
Numeric_2_Before=claw.frames[10].q[0,:]

## Plotting
fig, ax = plt.subplots()
ax.plot(xx,q_0_xi,'og',label=r"Semianalytic")
ax.plot(claw.frames[10].grid.p_centers[0],claw.frames[10].q[0,:],'-r', label=r"CLAWPACK")
ax.legend(loc='lower left') # upper left corner
ax.set_xlabel(r'Interval')
ax.set_ylabel(r'Values')
#ax.set_title('title');
plt.savefig("Comparision.jpg")

## Order of the convergence based on LeVeque's SIAM 2007 Book (Section A.5-A.6)
E_h_per_2=(h2*np.sum(np.abs(Numeric_2_Before-Semi_2_Before)**q))**(1.0/q)
R_h=(E_h)/(E_h_per_2)
p=math.log(R_h,2)
print p # order
