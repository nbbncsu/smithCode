def heat_fun_cu(q,a,b,L,k_cu,u_amb_cu):
   
   #
   # Function employed in the optimization routine for Example 12.1.
   #
   
   import math
   import numpy as np
   
   # global data_cu
   # global xdata
   #
   # Load the data and construct the x datapoints.
   #

   final_cu_data = np.loadtxt('final_cu_data.txt')

   data_cu = np.array(final_cu_data[1:16]);
   xdata =np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]); 
  
   h = q[0];
   Q = q[1];
   
   #
   #  Construct constants and solution
   #
   
   gamma_cu = math.sqrt(2*(a+b)*h/(a*b*k_cu));
   f1_cu = math.exp(gamma_cu*L)*(h + k_cu*gamma_cu);
   f2_cu = math.exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
   f3_cu = f1_cu/(f2_cu + f1_cu);
   c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
   c2_cu = Q/(k_cu*gamma_cu) + c1_cu;
  
   uvals_cu = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu; 
   res_cu = data_cu - uvals_cu;
   J = np.dot(res_cu,res_cu.T);

   return J




