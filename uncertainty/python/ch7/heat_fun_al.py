def  heat_fun_al(q,a,b,L,k_al,u_amb_al):

   import scipy as sp
   import numpy as np
   import math
 
   #   global data_al xdata
   final_al_data = np.loadtxt('final_al_data.txt')

   data_al = np.array(final_al_data[1:16]);
   xdata =np.array([10,14,18,22,26,30,34,38,42,46,50,54,58,62,66]);
 
   h = q[0]
   Q = q[1]
 
   #
   # Construct constants and solution
   # 
 
   gamma_al = math.sqrt(2*(a+b)*h/(a*b*k_al))
   f1_al = math.exp(gamma_al*L)*(h+k_al*gamma_al)
   f2_al = math.exp(-gamma_al*L)*(h - k_al*gamma_al)
   f3_al = f1_al/(f2_al + f1_al)
   c1_al = -Q*f3_al/(k_al*gamma_al)
   c2_al = Q/(k_al*gamma_al) + c1_al
 
   uvals_al = np.inner(c1_al,np.exp(-gamma_al*xdata)) + np.inner(c2_al,np.exp(gamma_al*xdata)) + u_amb_al
   res_al = np.array([data_al - uvals_al])
   J = np.dot(res_al,res_al.T)
   return J



