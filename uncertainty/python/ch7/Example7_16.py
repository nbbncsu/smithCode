#	Example7_16.py
#
# Code computes optimal parameters, the sensitivity matrices, and
# covariance matrix for Example 7.16.
#
# Required functions: heat_fun_al.m
# Required data: final_al_data.txt
#

#
# Import all required libraries
#

import scipy as sp
import scipy.optimize as opt
import numpy as np
import numpy.linalg as la
import math
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import heat_fun_al

global data_al
global xdata

#
# Load the data and construct the x datapoints.
#

final_al_data = np.loadtxt('final_al_data.txt')

data_al = np.array(final_al_data[1:16]);
xdata = np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
x0=1;xf=70;xstep=0.1; xvals = np.linspace(x0,xf,(xf-x0)/xstep+1);
u_amb = final_al_data[16];

#
# Input dimensions and material constants
# 

a = 0.95       #cm
b = 0.95       #cm
L = 70.0       #cm
k = 2.37       #W/cm C
n = 15         # Number of measurements
p = 2          # Number of parameters

#
# Optimize parameters and construct increments. The representations are 
# truncated to two significant digits to remain consistent with the 
# reported values.
#

h_init = 0.00183
Q_init = -15.93
q_init = np.array([h_init,Q_init])

modelfun = lambda q: heat_fun_al.heat_fun_al(q,a,b,L,k,u_amb)
q_opt,fval,iter,funcalls,warnflag = opt.fmin(modelfun,q_init,full_output=True)

h = q_opt[0]
Q = q_opt[1]

h = 0.00191
Q = -18.41

dh = 1e-10
dQ = 1e-4
h_p = h + dh
Q_p = Q + dQ

# 
# Construct constants and solution
# 

gamma = math.sqrt(2*(a+b)*h/(a*b*k))
gamma_h = (1/(2*h))*gamma
f1 = math.exp(gamma*L)*(h+k*gamma)
f2 = math.exp(-gamma*L)*(h-k*gamma)
f3 = f1/(f2+f1)
f1_h = math.exp(gamma*L)*(gamma_h*L*(h+k*gamma)+1+k*gamma_h)
f2_h = math.exp(-gamma*L)*(-gamma_h*L*(h-k*gamma) + 1 - k*gamma_h)
c1 = -Q*f3/(k*gamma)
c2 = Q/(k*gamma)+c1
f4 = Q/(k*gamma*gamma)
den2 = (f1+f2)**2
f3_h = (f1_h*(f1+f2) - f1*(f1_h+f2_h))/den2;
c1_h = f4*gamma_h*f3-(Q/(k*gamma))*f3_h
c2_h = -f4*gamma_h+c1_h
c1_Q = -(1/(k*gamma))*f3
c2_Q = (1/(k*gamma)) + c1_Q

gamma_hp = math.sqrt(2*(a+b)*h_p/(a*b*k))
f1_hp = math.exp(gamma_hp*L)*(h_p+k*gamma_hp)
f2_hp = math.exp(-gamma_hp*L)*(h_p - k*gamma_hp)
f3_hp = f1_hp/(f2_hp + f1_hp)
c1_hp = -Q*f3_hp/(k*gamma_hp)
c2_hp = Q/(k*gamma_hp) + c1_hp

c1_Qp = -Q_p*f3/(k*gamma)
c2_Qp = Q_p/(k*gamma) + c1_Qp

uvals = c1*np.exp(-gamma*xvals) + c2*np.exp(gamma*xvals) + u_amb
uvals_data = c1*np.exp(-gamma*xdata) + c2*np.exp(gamma*xdata) + u_amb
uvals_Q_data = c1_Q*np.exp(-gamma*xdata) + c2_Q*np.exp(gamma*xdata)
uvals_h_data = c1_h*np.exp(-gamma*xdata) + c2_h*np.exp(gamma*xdata) + gamma_h*np.multiply(xdata,(-c1*np.exp(-gamma*xdata) + c2*np.exp(gamma*xdata)))

uvals_data_hp = c1_hp*np.exp(-gamma_hp*xdata) + c2_hp*np.exp(gamma_hp*xdata) +u_amb
uvals_data_Qp = c1_Qp*np.exp(-gamma*xdata) + c2_Qp*np.exp(gamma*xdata) + u_amb

S_d = np.concatenate((np.array([np.ones(15)]), np.array([xdata]),np.array([xdata**2]), np.array([xdata**3]), np.array([xdata**4])),axis=0)

h_fd = np.array([np.multiply((1/dh),(uvals_data_hp - uvals_data))])
Q_fd = np.array([np.multiply((1/dQ),(uvals_data_Qp - uvals_data))])

res = np.array([data_al-uvals_data])

#
# Construct the analytic and finite difference sensitivity matrices
#

sens_mat = np.concatenate((np.array([uvals_Q_data]), np.array([uvals_h_data])),axis=0)
sens_mat_fd = np.concatenate((Q_fd, h_fd),axis=0)

#
# Construct the measurement covariance sigma2 and the covariance matrices
# V and V_fd constructed using the analytic and finite difference 
# sensitivity relations.
#

sigma2 = np.multiply((1/(n-p)),np.dot(res,res.T))
V = sigma2*la.inv(np.dot(sens_mat_fd, sens_mat_fd.T))

Sens_mat_aug = np.concatenate((sens_mat, S_d),axis=0)
Umat,Smat,Vmat = la.svd(Sens_mat_aug,full_matrices=1,compute_uv=1)

#
# Construct the 95% confidence intervals.
#

tval = 2.1604
int_Q = np.array([Q-math.sqrt(V[0,0])*tval,Q + math.sqrt(V[0,0])*tval])
int_h = np.array([h - math.sqrt(V[1,1])*tval, h + math.sqrt(V[1,1])*tval])

# 
# Plot solutions. In Figures 3 and 4, we plot the sensitivites obtained both
# analytically and using finite differences to show that to within visual 
# accuracy, they are the same.
#

plt.figure(1)
plt.plot(xvals,uvals,'k',linewidth=2.0,label='Model')
plt.plot(xdata,data_al,'o',linewidth=5.0,label='Data')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (^oC)')
plt.legend()

plt.figure(2)
plt.plot(xdata,res.T,'o',linewidth=6.0)
plt.plot(xvals,np.multiply(0,xvals),linewidth=2.0)
plt.xlabel('Distance (cm)')
plt.ylabel('Residuals (^oC)')

plt.figure(3)
plt.plot(xdata,uvals_Q_data,'x',ms=10)
plt.plot(xdata,Q_fd.T,'o',ms=6)
plt.xlabel('Distance (cm)')
plt.ylabel('Sensitivity in \Phi')

plt.figure(4)
plt.plot(xdata,uvals_h_data,'x',ms=10)
plt.plot(xdata,h_fd.T,'o', ms=6)
plt.xlabel('Distance (cm)')
plt.ylabel('Sensitivity in h')
plt.show()







