#
#                       heat_code_mcmc.py
#
#
# This code illustrates the implementation of the Metropolis algorithm 8.5
# used to construct the chains in Figure 8.10 for the heat Example 8.12.
# The code takes approximately 30 seconds to run.  This can be reduced
# by reducing the number of chain iterations.
#

#
# Required functions: kde.m
# Required data: final_al_data.txt
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
# Load data
#

final_al_data = np.loadtxt('final_al_data.txt')

data = final_al_data[1:16];
xdata = np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
x0 = 10; x1 = 70; xstep = 0.1; xnum=(x1-x0)/xstep+1
xvals = np.linspace(x0,x1,xnum);
u_amb = final_al_data[16];

#
# Input dimensions and material constants
#

a = 0.95;   # cm
b = 0.95;   # cm
L = 70.0;   # cm
k = 2.37;   # W/cm C
h = 0.00191;
Q = -18.41; 
n = float(15);
p = 2;

#
# Construct constants and solution
#

gamma = math.sqrt(2*(a+b)*h/(a*b*k));
gamma_h = (1/(2*h))*gamma;
f1 = math.exp(gamma*L)*(h + k*gamma);
f2 = math.exp(-gamma*L)*(h - k*gamma);
f3 = f1/(f2 + f1);
f1_h = math.exp(gamma*L)*(gamma_h*L*(h+k*gamma) + 1 + k*gamma_h);
f2_h = math.exp(-gamma*L)*(-gamma_h*L*(h-k*gamma) + 1 - k*gamma_h);
c1 = -Q*f3/(k*gamma);
c2 = Q/(k*gamma) + c1;
f4 = Q/(k*gamma*gamma);
den2 = (f1+f2)**2;
f3_h = (f1_h*(f1+f2) - f1*(f1_h+f2_h))/den2;
c1_h = f4*gamma_h*f3 - (Q/(k*gamma))*f3_h;
c2_h = -f4*gamma_h + c1_h;
c1_Q = -(1/(k*gamma))*f3;
c2_Q = (1/(k*gamma)) + c1_Q;

uvals = c1*np.exp(-gamma*xvals) + c2*np.exp(gamma*xvals) + u_amb;
uvals_data = c1*np.exp(-gamma*xdata) + c2*np.exp(gamma*xdata) + u_amb;
uvals_Q_data = c1_Q*np.exp(-gamma*xdata) + c2_Q*np.exp(gamma*xdata);
uvals_h_data = c1_h*np.exp(-gamma*xdata) + c2_h*np.exp(gamma*xdata) + gamma_h*np.multiply(xdata,(-c1*np.exp(-gamma*xdata) + c2*np.exp(gamma*xdata)));


res = data - uvals_data;

sens_mat = np.concatenate(([uvals_Q_data],[uvals_h_data]), axis=0);
sigma2 = (1/(n-p))*(np.inner(res,res.T))
V = sigma2*la.inv(np.dot(sens_mat,sens_mat.T));

#
# Set MCMC parameters
#

N = 10000;
R = la.cholesky(V);
q_old = np.array([[Q],[h]]);
SS_old = np.inner(res,res.T);
n0 = 0.001;
sigma02 = sigma2;
aval = 0.5*(n0 + 15);
bval = 0.5*(n0*sigma02 + SS_old);
sigma2 = 1/np.random.gamma(aval,1/bval);
accept = 0;

#
#  Run the Metropolis algorithm for N iterations.
#

Q_MCMC = np.zeros(shape=(2,int(N)))
Sigma2 = np.zeros(shape=(int(N),1))
  
for i in range(0,int(N)):
   z = np.random.randn(2,1); 
   q_new = q_old + np.dot(R,z);
   Q = q_new[0,0];
   h = q_new[1,0];
   gamma = math.sqrt(2*(a+b)*h/(a*b*k));
   f1 = math.exp(gamma*L)*(h + k*gamma);
   f2 = math.exp(-gamma*L)*(h - k*gamma);
   f3 = f1/(f2 + f1);
   c1 = -Q*f3/(k*gamma);
   c2 = Q/(k*gamma) + c1;
   uvals_data = c1*np.exp(-gamma*xdata) + c2*np.exp(gamma*xdata) + u_amb;
   res = data - uvals_data;
   SS_new = np.inner(res,res.T);
   u_alpha = np.random.random(1);
   term = np.exp(-.5*(SS_new-SS_old)/sigma2);
   alpha = min(1,term);
   if u_alpha < alpha:
      Q_MCMC[0,i] = Q
      Q_MCMC[1,i] = h
      q_old = q_new;
      SS_old = SS_new;
      accept = accept + 1;
   else:
      Q_MCMC[0,i] = q_old[0,0]
      Q_MCMC[1,i] = q_old[1,0]
   Sigma2[i] = sigma2;
   bval = 0.5*(n0*sigma02 + SS_old);
   sigma2 = 1/np.random.gamma(aval,1/bval);

Qvals = Q_MCMC[0,0:];
hvals = Q_MCMC[1,0:];
  
#
# Use kde to construct densities for Q and h.
#

range_Q = max(Qvals) - min(Qvals);
range_h = max(hvals) - min(hvals);
Q_min = min(Qvals)-range_Q/10;
Q_max = max(Qvals)+range_Q/10;
h_min = min(hvals)-range_h/10;
h_max = max(hvals)+range_h/10;
density1 = stats.gaussian_kde(Qvals)
Qmesh = np.linspace(Q_min,Q_max,16384)
density_Q = density1.evaluate(Qmesh)
density2 = stats.gaussian_kde(hvals)
hmesh = np.linspace(h_min,h_max,16384)
density_h = density2.evaluate(hmesh)


print accept/N

#
# Plot solutions
#
 
plt.figure(1)
plt.plot(Qvals,'-')
plt.axis([0,N,-19.3,-17.5])
plt.xlabel('Chain Iteration')
plt.ylabel('Parameter Q')

plt.figure(2)
plt.plot(hvals,'-')
plt.axis([0,N,1.84e-3,2e-3])
plt.xlabel('Chain Iteration')
plt.ylabel('Parameter h')

plt.figure(3)
plt.plot(Sigma2)
plt.title('Measurement Error Variance \sigma^2')

plt.figure(4)
plt.plot(Qmesh,density_Q,'k-')
plt.axis([-19.5,-17.5,0,3])
plt.xlabel('Parameter Q')

plt.figure(5)
plt.plot(hmesh,density_h,'k-')
plt.axis([1.8e-3,2e-3,0,3e4])
plt.xlabel('Parameter h')

plt.figure(6)
plt.scatter(Qvals,hvals)
plt.axis([-19.2,-17.6,1.83e-3,2e-3])
plt.xlabel('Parameter Q')
plt.ylabel('Parameter h')

plt.show()
