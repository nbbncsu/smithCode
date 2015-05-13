def heat_fun_cu_quad(q,a,b,L,k_cu,u_amb_cu):

   #
   # Function employed in the optimization routine for Example 12.6.
   #
   # global data_cu xdata
   import math
   import numpy as np

   final_cu_data = np.loadtxt('final_cu_data.txt')
   data_cu = final_cu_data[1:16]
   data_cu = np.delete(data_cu,10)

   xdata = np.array([10,14,18,22,26,30,34,38,42,46,54,58,62,66]);

   h = q[0];
   Q = q[1];
   b0 = q[2];
   b1 = q[3];
   b2 = q[4];

   #
   #  Construct constants and solution
   #

   gamma_cu = np.sqrt(2*(a+b)*h/(a*b*k_cu));
   f1_cu = math.exp(gamma_cu*L)*(h + k_cu*gamma_cu);
   f2_cu = math.exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
   f3_cu = f1_cu/(f2_cu + f1_cu);
   c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
   c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

   uvals_cu = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu + b0 + b1*xdata + b2*(xdata**2); 
   res_cu = np.array([data_cu - uvals_cu]);
   J = np.dot(res_cu,res_cu.T);

   return J



