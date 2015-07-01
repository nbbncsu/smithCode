#
#                       spring_mcmc_C_K_sigma.py
#
#
# This code illustrates the implementation of the Metropolis algorithm
# in Example 8.7, Case ii.  Here we are sampling the stiffness and damping
# parameter q = [C,K] as well as the variance sigma^2 for the measurement
# error. 
#

#
# Import all required libraries
#

import math
import scipy as sp
import numpy as np
import numpy.linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import bigfloat

#
# Input true stiffness and damping parameters and error variance
#

K = 20.5;
C = 1.5;
sigma = .1;
var = sigma**2;
N = 10000;

factor = math.sqrt(K - C**2/4);
denom = math.sqrt(4*K - C**2);

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
dydc = np.array([np.multiply(np.exp(-C*t/2),(np.multiply((C*t/denom),np.sin(factor*t)) - np.multiply(t,np.cos(factor*t))))])

dydk = np.array([np.multiply((-2*t/denom),np.multiply(np.exp(-C*t/2),np.sin(factor*t)))])

error = sigma*np.random.randn(n)
obs = y + error

#
# Construct the sensitivity matrix chi and covariance matrix V.
#

chi = np.concatenate((dydc.T, dydk.T), axis=1)
p = np.dot(chi.T,chi)
V = la.inv(p)*var;
R = la.cholesky(V);

q_old = np.array([[C],[K]])
SS_old = np.inner((obs-y),(obs-y).T);

#
# Construct the parameters used when sampling the error variance in the Metropolis algorithm.
#

n0 = .001;
sigma02 = var;
aval = 0.5*(n0 + n);
bval = 0.5*(n0*sigma02 + SS_old);
sigma2 = 1/np.random.gamma(aval,1/bval);

#
# Construct a Metropolis chain of length N
#

Q_MCMC = np.zeros(shape=(2,int(N)))
Sigma2 = np.zeros(shape=(int(N),1))

for i in range(0,int(N)):
   z = np.random.randn(2,1); 
   q_new = q_old + np.dot(R,z);
   C_new = q_new[0,0];
   K_new = q_new[1,0];
   factor = math.sqrt(K_new - C_new**2/4);
   y_new = np.multiply(2*np.exp(-C_new*t/2),np.cos(factor*t))
   u_alpha = np.random.random(1);
   SS_new = np.inner((obs-y_new),(obs-y_new).T);
   term = np.exp(-.5*(SS_new-SS_old)/sigma2);
   alpha = min(1,term);
   if u_alpha < alpha:
      Q_MCMC[0,i] = C_new
      Q_MCMC[1,i] = K_new
      q_old = q_new;
      SS_old = SS_new;
   else:
      Q_MCMC[0,i] = q_old[0,0];
      Q_MCMC[1,i] = q_old[1,0];
   Sigma2[i] = sigma2;
   bval = 0.5*(n0*sigma02 + SS_old);
   sigma2 = 1/np.random.gamma(aval,1/bval);

Cvals = Q_MCMC[0,0:];
Kvals = Q_MCMC[1,0:];

#
# Use kde to construct densities for C and K
#

range_C = max(Cvals) - min(Cvals);
range_K = max(Kvals) - min(Kvals);
C_min = min(Cvals)-range_C/10;
C_max = max(Cvals)+range_C/10;
K_min = min(Kvals)-range_K/10;
K_max = max(Kvals)+range_K/10;
density = stats.gaussian_kde(Cvals.T);
Cmesh = np.linspace(C_min,C_max,16384)
density_C = density.evaluate(Cmesh)
# bandwidth_C, Cmesh, cdf_C    #16384 is hardcoded from before
density2 = stats.gaussian_kde(Kvals);
Kmesh = np.linspace(K_min,K_max,16384)
density_K = density2.evaluate(Kmesh)
# bandwidth_K, Kmesh, cdf_K

#
# Plot results.  Note that the axis limits may change as a function of
# each new generated data set.
#

plt.figure(1)
plt.plot(Cvals,'-')
plt.axis([0,N,1.43,1.57])
plt.xlabel('Chain Iteration')
plt.ylabel('Damping Parameter C')

plt.figure(2)
plt.plot(Kvals,'-')
#  axis([0,N,20.2,20.9])
plt.xlabel('Chain Iteration')
plt.ylabel('Stiffness Parameter K')

plt.figure(3)
plt.plot(Sigma2)
plt.title('Measurement Error Variance \sigma^2')

plt.figure(4)
plt.plot(Cmesh,density_C,'k-')
plt.axis([1.4,1.6,0,25])
plt.xlabel('Damping Parameter C')

plt.figure(5)
plt.plot(Kmesh,density_K,'k-')
#  axis([20,21,0,5])
plt.xlabel('Stiffness Parameter K')

plt.figure(6)
plt.scatter(Cvals,Kvals)
# box on
# axis([-19.2,-17.6,1.83e-3,2e-3])
plt.xlabel('Damping Parameter C')
plt.ylabel('Stiffness Parameter K')

plt.show()


