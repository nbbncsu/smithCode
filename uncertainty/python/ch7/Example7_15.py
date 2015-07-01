#	Example7_15.py
#
# Construct synthetic data and sampling distribution for the damping
# parameter C as illustrated in Example 7.15.
#
# Required function: spring_function.py
#

#
# Import all required libraries
#

import scipy as sp
import numpy as np
import math
import numpy.linalg as la
import scipy.stats as stats
import scipy.optimize as opt
import matplotlib as mpl
import matplotlib.pyplot as plt
import spring_function

#
# Specify true parameter and measurement variance values
#

m = 1
k = 20.5
c = 1.5
sigma = 0.1
var = sigma**2

#
# Construct synthetic data, sensitivity value, and covariance value V.
#

num = math.sqrt(4*m*k-c**2)
den = 2*m

t0=0; tf=5; tstep=0.01;n=(tf-t0)/tstep+1;
t=np.linspace(t0,tf,n)
n = len(t)
error = sigma*np.random.randn(n)

y = np.multiply(2*np.exp((-c/den)*t),np.cos((num/den)*t))
dydc = np.array([np.multiply(-(t/m),np.multiply(np.exp((-c/den)*t),np.cos((num/den)*t)))+np.multiply((c*t/(m*num)),np.multiply(np.exp((-c/den)*t),np.sin((num/den)*t)))])
print dydc.shape
print dydc.T.shape
V = la.inv(np.dot(dydc, np.transpose(dydc)))*var;

obs = y+error

#
# Construct sampling distribution and determine optimal value of c
# for this data.

C0=1.4; Cf=1.6; Cstep=0.001;N=(Cf-C0)/Cstep+1;
C=np.linspace(C0,Cf,N)
sample_dist=stats.norm.pdf(C,1.5,np.sqrt(V))
print C.shape
print sample_dist.shape

modelfun = lambda c: spring_function.spring_function(c,m,k,obs,t)
c_opt,fval,iter,funcalls,warnflag=opt.fmin(modelfun,c,full_output=True)

#
# Plot results
#
plt.figure(1)
plt.plot(t,y,'k',label='Model')
plt.plot(t,obs,'x',label='Synthetic Data')
plt.plot(t,0*y,'k')
plt.xlabel('Time (s)')
plt.ylabel('Displacement')
plt.legend()

plt.figure(2)
plt.plot(t,error,'x')
plt.plot(t,0*y,'k')
plt.plot(t,2*sigma*np.ones(len(t)),'--r')
plt.plot(t,-2*sigma*np.ones(len(t)),'--r')
plt.xlabel('Time (s)')
plt.ylabel('Residuals')

plt.figure(3)
plt.plot(C,sample_dist.T)
plt.xlabel('Damping Parameter c')
plt.ylabel('Sampling Density')
plt.show()






