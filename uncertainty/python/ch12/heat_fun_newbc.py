def heat_fun_newbc(q,a,b,L,k_al,k_cu,u_amb_cu,u_amb_al):
 
   #global data_cu data_al xdata
   import math
   import numpy as np

   final_cu_data = np.loadtxt('final_cu_data.txt')
   final_al_data = np.loadtxt('final_al_data.txt')

   data_cu = final_cu_data[1:16];
   data_al = final_al_data[1:16];
   xdata = np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
   u_amb_cu = final_cu_data[16];
   u_amb_al = final_al_data[16];
   
   h = q[0];
   Tsource = q[1];
   eta = q[2];

   #
   #  Construct constants and solution
   #
 
   gamma_cu = np.sqrt(2*(a+b)*h/(a*b*k_cu));
   Gamma_cu = (k_cu*gamma_cu + eta)/(k_cu*gamma_cu - eta);
   f1_cu = math.exp(gamma_cu*L)*Gamma_cu*(h + k_cu*gamma_cu);
   f2_cu = math.exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
   f3_cu = math.exp(gamma_cu*L)*(h + k_cu*gamma_cu);
   f4_cu = f3_cu/(f2_cu + f1_cu);
   c1_cu = f4_cu*eta*(u_amb_cu - Tsource)/(eta - k_cu*gamma_cu);
   c2_cu = eta*(Tsource - u_amb_cu)/(eta - k_cu*gamma_cu) + c1_cu*Gamma_cu;
 
   gamma_al = np.sqrt(2*(a+b)*h/(a*b*k_al));
   Gamma_al = (k_al*gamma_al + eta)/(k_al*gamma_al - eta);
   f1_al = math.exp(gamma_al*L)*Gamma_al*(h + k_al*gamma_al);
   f2_al = math.exp(-gamma_al*L)*(h - k_al*gamma_al);
   f3_al = math.exp(gamma_al*L)*(h + k_al*gamma_al);
   f4_al = f3_al/(f2_al + f1_al);
   c1_al = f4_al*eta*(u_amb_al - Tsource)/(eta - k_al*gamma_al);
   c2_al = eta*(Tsource - u_amb_al)/(eta - k_al*gamma_al) + c1_al*Gamma_al;
 
   uvals_cu = c1_cu*np.exp(-gamma_cu*xdata) + c2_cu*np.exp(gamma_cu*xdata) + u_amb_cu;
   uvals_al = c1_al*np.exp(-gamma_al*xdata) + c2_al*np.exp(gamma_al*xdata) + u_amb_al;
 
   res_cu = np.array([data_cu - uvals_cu]);
   res_al = np.array([data_al - uvals_al]);
   J = np.dot(res_cu,res_cu.T) + np.dot(res_al,res_al.T);


   return J
