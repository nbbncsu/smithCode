#
#                           Example12_2.m
#
#
# Code computes and plots the optimal parameters for the heat model (3.21)
# simultaneously using copper and aluminum data. 
#
# Required functions: heat_fun.m
# Required data: final_cu_data.txt
#                final_al_data.txt


#global data_cu data_al xdata 

import math
import scipy as sp
import scipy.optimize as opt
import numpy as np
import numpy.linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import heat_fun
  
#
# Load the data and construct the x datapoints.
#

final_cu_data = np.loadtxt('final_cu_data.txt')
final_al_data = np.loadtxt('final_al_data.txt')

data_cu = final_cu_data[1:16];
data_al = final_al_data[1:16];
xdata = np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
u_amb_cu = final_cu_data[16];
u_amb_al = final_al_data[16];
x0=0;xf=70;xstep=0.1; xvals = np.linspace(x0,xf,(xf-x0)/xstep+1);


#
# Input dimensions and material constants
#

a = 0.95;      # cm
b = 0.95;      # cm
L = 70.0;      # cm
k_cu = 4.01;   # W/cm C
k_al = 2.37;   # W/cm C
n = 15;        # Number of measurements
p = 2;         # Number of parameters
  
h_init = 0.00183;
Q_init = -15.93;
q_init = np.array([h_init,Q_init]);

#
# Optimize parameters
#

modelfun = lambda q: heat_fun.heat_fun(q,a,b,L,k_al,k_cu,u_amb_al,u_amb_cu);
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

gamma_al = math.sqrt(2*(a+b)*h/(a*b*k_al));
f1_al = math.exp(gamma_al*L)*(h + k_al*gamma_al);
f2_al = math.exp(-gamma_al*L)*(h - k_al*gamma_al);
f3_al = f1_al/(f2_al + f1_al);
c1_al = -Q*f3_al/(k_al*gamma_al);
c2_al = Q/(k_al*gamma_al) + c1_al;

uvals_cu = c1_cu*np.exp(-gamma_cu*xvals) + c2_cu*np.exp(gamma_cu*xvals) + u_amb_cu;
uvals_al = c1_al*np.exp(-gamma_al*xvals) + c2_al*np.exp(gamma_al*xvals) + u_amb_al;

uvals_cu_data = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu;
uvals_al_data = c1_al*np.exp(-gamma_al*xdata) + c2_al*np.exp(gamma_al*xdata) + u_amb_al;

res_cu = np.array([data_cu - uvals_cu_data]);
res_al = np.array([data_al - uvals_al_data]);
  
#
# Plot the results
#

plt.figure(1)
plt.plot(xvals,uvals_cu,linewidth=2.0,label='Model')
plt.axis([0,70,20,uvals_cu[0]])
plt.plot(xdata,data_cu,'o',linewidth=5.0,label='Data')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (^oC)')
plt.legend()

plt.figure(2)
plt.plot(xvals,uvals_al,linewidth=2.0,label='Model')
plt.axis([0,70,20,uvals_al[0]])
plt.plot(xdata,data_al,'o',linewidth=5.0,label='Data')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (^oC)')
plt.legend()

plt.show()
  
 
  
