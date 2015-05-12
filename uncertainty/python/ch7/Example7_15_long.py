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
import scipy.stats as stats
import scipy.optimize as opt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy.linalg as la
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
dydc = np.array([np.multiply(-(t/m),np.multiply(np.exp((-c/den)*t),np.cos((num/den)*t)))+np.multiply((c*t/(m*num)),np.multiply(np.exp((-c/den)*t),np.sin((num/den)*t)))]);
V = la.inv(np.dot(dydc, dydc.T))*var;
obs = y+error

#
# Construct sampling distribution and determine optimal value of c
# for this data.

C0=1.4; Cf=1.6; Cstep=0.001;N=(Cf-C0)/Cstep+1;
C=np.linspace(C0,Cf,N)
sample_dist=stats.norm.pdf(C,1.5,np.sqrt(V))

modelfun = lambda c: spring_function.spring_function(c,m,k,obs,t)
c_opt,fval,iter,funcalls,warnflag=opt.fmin(modelfun,c,full_output=True)
conf_int = np.array([c_opt-1.96*np.sqrt(V), c_opt+1.96*np.sqrt(V)])

#
# Construct the density by repeated sampling. The final KDE is performed
# using ksdensity. The variable count compiles the number of 95% confidence
# intervals that contain the true parameter value.
#

Cvals = np.zeros(10000)

count = 0
for j in range(1,10000):
	error = sigma*np.random.randn(len(t))
	y = np.multiply(2*np.exp((-c/den)*t),np.cos((num/den)*t))
	obs = y + error
	modelfun = lambda c: spring_function.spring_function(c,m,k,obs,t)
	c_opt,fval,iter,funcalls,warnflag=opt.fmin(modelfun,c,full_output=True,disp=False)
	c_lower = c_opt - 1.96*np.sqrt(V)
	c_upper = c_opt + 1.96*np.sqrt(V)
	if (c_lower <= c) and (c<= c_upper):
		count = count+1
	Cvals[j] = c_opt
count

#[density_C, Cmesh] = stats.gaussian_kde(Cvals)
density = stats.gaussian_kde(Cvals.T)
Cmesh = np.linspace(Cvals.min(), Cvals.max(), 1000)
density_C = density.evaluate(Cmesh)

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
plt.axis([1.4,1.6,0,25])
plt.xlabel('Damping Parameter c')
plt.ylabel('Sampling Density')

plt.figure(4)
plt.plot(Cmesh,density_C,'k--',label='Constructed Density')
plt.plot(C,sample_dist.T,label='Sampling Density')
plt.axis([1.4,1.6,0,25])
plt.xlabel('Optimal C')
plt.ylabel('Density')
plt.legend()
plt.show()





