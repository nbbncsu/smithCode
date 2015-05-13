#
#                           Example12_6.m
#
#
# Code computes and plots the optimal parameters for the heat model (3.21)
# using the quadratic discrepancy relation (12.5).  The 11th data point is
# deleted as an outlier.
#
# Required functions: heat_fun_cu_quad.m
# Required data: final_cu_data.txt
#

# global data_cu xdata

import math
import scipy as sp
import scipy.optimize as opt
import numpy as np
import numpy.linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import heat_fun_cu_quad
  
#
# Load the data and construct the x datapoints.
#

final_cu_data = np.loadtxt('final_cu_data.txt')
data_cu = final_cu_data[1:16];
data_cu = np.delete(data_cu,10)
xdata = np.array([10,14,18,22,26,30,34,38,42,46,54,58,62,66]);
u_amb_cu = final_cu_data[16];
x0=0;xf=70;xstep=0.1; xvals = np.linspace(x0,xf,(xf-x0)/xstep+1);

#
# Input dimensions and material constants
#

a = 0.95;      # cm
b = 0.95;      # cm
L = 70.0;      # cm
k_cu = 4.01;   # W/cm C
n = 14;        # Number of measurements
p = 2;         # Number of parameters

h_init = 0.00183;
Q_init = -15.93;
b0_init = .0;
b1_init = 0.0;
b2_init = 0.0;
q_init = [h_init,Q_init,b0_init,b1_init,b2_init];

#
# Optimize parameters
#

modelfun = lambda q: heat_fun_cu_quad.heat_fun_cu_quad(q,a,b,L,k_cu,u_amb_cu);
q_opt,fval,iter,funcalls,warnflag = opt.fmin(modelfun,q_init,full_output=True,disp=False);

print q_opt 


h = q_opt[0];
Q = q_opt[1];
b0 = q_opt[2];
b1 = q_opt[3];
b2 = q_opt[4];
  
#
# Compute solution using optimal parameter values.
#

gamma_cu = np.sqrt(2*(a+b)*h/(a*b*k_cu));
f1_cu = math.exp(gamma_cu*L)*(h + k_cu*gamma_cu);
f2_cu = math.exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
f3_cu = f1_cu/(f2_cu + f1_cu);
c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

uvals_cu = c1_cu*np.exp(-gamma_cu*xvals) + c2_cu*np.exp(gamma_cu*xvals) + u_amb_cu + b0 + b1*xvals + b2*(xvals**2);
uvals_cu_data = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu + b0 + b1*xdata + b2*(xdata**2);

res_cu = np.array([data_cu - uvals_cu_data]);
  
#
# Plot the results
#

plt.figure(1)
plt.plot(xvals,uvals_cu,linewidth=2.0,label='Model')
plt.axis([0,70,20,uvals_cu[1]])
plt.plot(xdata,data_cu,'o',linewidth=5.0,label='Data')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (^oC)')
plt.legend()


plt.figure(2)
plt.plot(xdata,res_cu.T,'o',linewidth=6.0)
plt.axis([0,70,-.12,.2])
plt.plot(xvals,0*xvals,linewidth=2.0)
plt.xlabel('Distance (cm)')
plt.ylabel('Residuals (^oC)')

plt.show()  


