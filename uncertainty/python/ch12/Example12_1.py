#
#                           Example12_1.m
#
#
# Code computes and plots the optimal parameters for the heat model (3.21)
# using copper data.
#
# Required functions: heat_fun_cu.m
# Required data: final_cu_data.txt
#


import math
import scipy as sp
import scipy.optimize as opt
import numpy as np
import numpy.linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import heat_fun_cu


global data_cu
global xdata
  
#
# Load the data and construct the x datapoints.
#
  
final_cu_data = np.loadtxt('final_cu_data.txt')

data_cu = np.array(final_cu_data[1:16]);
xdata = np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
u_amb_cu = final_cu_data[16];
x0=0;xf=70;xstep=0.1; xvals = np.linspace(x0,xf,(xf-x0)/xstep+1);

#
# Input dimensions and material constants
#

a = 0.95;      # cm
b = 0.95;      # cm
L = 70.0;      # cm
k_cu = 4.01;   # W/cm C
n = 15;        # Number of measurements
p = 2;         # Number of parameters

h_init = 0.00183;
Q_init = -15.93;
q_init = np.array([h_init,Q_init]);

#
# Optimize parameters
#

modelfun = lambda q: heat_fun_cu.heat_fun_cu(q,a,b,L,k_cu,u_amb_cu);
q_opt,fval,iter,funcalls,warnflag = opt.fmin(modelfun,q_init,full_output=True);

print q_opt 

h = q_opt[0];
Q = q_opt[1];

#
# Compute solution using optimal parameter values.
#

gamma_cu = math.sqrt(2*(a+b)*h/(a*b*k_cu));
f1_cu = math.exp(gamma_cu*L)*(h + k_cu*gamma_cu);
f2_cu = math.exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
f3_cu = f1_cu/(f2_cu + f1_cu);
c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

uvals_cu = c1_cu*np.exp(-gamma_cu*xvals) + c2_cu*np.exp(gamma_cu*xvals) + u_amb_cu;
uvals_cu_data = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu;

res_cu = data_cu - uvals_cu_data;

#
# Plot the results
#

plt.figure(1)
plt.plot(xvals,uvals_cu,label='Model')
plt.axis([0,70,20,uvals_cu[1]])
plt.plot(xdata,data_cu,'o',label='Data')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (^oC)')
plt.legend()

plt.figure(2)
plt.plot(xdata,res_cu,'o',linewidth=6.0)
plt.axis([0,70,-.4,.45])
plt.plot(xvals,0*xvals)
plt.xlabel('Distance (cm)')
plt.ylabel('Residuals (^oC)')
plt.show()  


