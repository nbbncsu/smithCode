# 	spring_mcmc_C.py
#
#
# This code illustrates the implementation of the Metropolis algorithm
# in Example 8.7, Case i.
#


# Import all required libraries
import math
import scipy as sp
import numpy as np
import numpy.linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import bigfloat

#
#Initialize all constants
#

K=20.5;
C=1.5;
C_old = C
C0 = C
sigma=0.1;
var=sigma**2;
N = 10000

factor = math.sqrt(K-C**2/4)
denom = math.sqrt(4*K-C**2)

#
# Construct solution on [0,5], the analytic sensitivity relations, and 
# the observations, which include the error.
# 

t0=0;
tf=5;
tstep=0.01;
n=(tf-t0)/tstep+1
t=np.linspace(t0,tf,n)

y = (np.multiply(2*np.exp(-C*t/2),np.cos(factor*t)))
y0 = y
dydc=np.multiply(np.exp(-C*t/2),(np.multiply(C*t/denom,np.sin(factor*t))-np.multiply(t,np.cos(factor*t))))

error=sigma*np.random.randn(n)
obs = y + error

#
# Construct the sensitivity matrix chi and covariance matrix V.
# 

chi=np.transpose(dydc)
V=(1/(np.inner(dydc,chi)))*var
R=math.sqrt(V)

SS_old = np.inner(obs-y,np.transpose(obs-y))

#
# Construct a Metropolis chain of N
# 

Cvals = np.zeros(shape=(int(N),1))

for i in range(0,int(N)):
   z = np.random.randn(1)
   C_new = C_old + R*z
   factor = math.sqrt(K-C_new**2/4)
   y_new = np.multiply(2*np.exp(np.multiply(-C_new/2,t)),np.cos(np.multiply(factor,t)))
   u_alpha = np.random.random(1)
   SS_new = np.inner(obs-y_new,np.transpose(obs-y_new))
   term = np.exp(-0.5*(SS_new-SS_old)/var)
   alpha = min(1,term)
   if u_alpha < alpha:
      Cvals[i] = C_new
      C_old = C_new
      SS_old = SS_new
   else:
      Cvals[i] = C_old

#
# Use kde to construct the density for C from the Metropolis chain.
# 

range_C = max(Cvals) - min(Cvals)
C_min = min(Cvals) - range_C/10
C_max = max(Cvals) + range_C/10
density = stats.gaussian_kde(Cvals.T)
C_mesh = np.linspace(C_min,C_max,16384)
density_C = density.evaluate(C_mesh)

#
# Construct the posterior density directly by approximating Bayes' rule using a
# a midpoint quadrature rule.
#

h=0.001; Nq=(0.2/h)
quadpt=np.zeros(shape=(int(Nq),1))
posterior=np.zeros(shape=(int(Nq),1))
C_post_axis=np.zeros(shape=(int(Nq),1))

for i in range(0,int(Nq)):
   C = 1.4 + (i-1)*h
   factor = math.sqrt(K-C**2/4)
   y = 2*np.multiply(np.exp(np.multiply(-C/2,t)),np.cos(np.multiply(factor,t)))
   SS=np.inner(obs-y,np.transpose(obs-y))
   for j in range(1,int(Nq)):
      Cq = 1.4 + (j-1)*h
      factor = math.sqrt(K-Cq**2/4)
      yq = 2*np.multiply(np.exp(np.multiply(-Cq/2,t)),np.cos(np.multiply(factor,t)))
      SSq = np.inner(obs-yq,np.transpose(obs-yq))
      quadpt[j] = np.exp(-0.5*(SSq-SS)/var)
   posterior[i]=1/(h*np.sum(quadpt))
   C_post_axis[i] = C

#
# Construct the asymptotic normal sampling distribution
# 

mean=h*np.sum(np.multiply(C_post_axis,posterior))
sample_dist=stats.norm.pdf(C_post_axis,C0,np.sqrt(V))

#
#Plot results.
#

plt.figure()
plt.plot(t,y0,linewidth=2.0)
plt.plot(t,np.multiply(0,y0),'k',linewidth=2.0)
plt.plot(t,obs,'x',linewidth=1.0)
plt.ylabel('Displacement')
plt.xlabel('Time (s)')
#plt.show()

plt.figure()
plt.plot(t,error,'x',linewidth=2.0)
plt.plot(t,np.multiply(0,y),'k',linewidth=2.0)
plt.plot(t,2*sigma*np.ones(shape=(len(t),1)),'--r',linewidth=2.0)
plt.plot(t,-2*sigma*np.ones(shape=(len(t),1)),'--r',linewidth=2.0)
plt.xlabel('Time (s)')
plt.ylabel('Residuals')

plt.figure()
plt.plot(Cvals,'-',linewidth=1.0)
plt.ylabel('Damping Parameter C')
plt.xlabel('Chain Iteration')
#plt.show()

plt.figure()
plt.plot(C_post_axis,posterior,'b-',linewidth=3.0,label='Bayes')
plt.plot(C_mesh,density_C,'g-.', linewidth=3.0,label='MCMC')
plt.plot(C_post_axis,sample_dist,'r--',linewidth=3.0,label='OLS')
plt.legend()
plt.xlabel('Damping Parameter C')
plt.show()



