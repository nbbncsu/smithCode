def spring_function(c,m,k,obs,t):
   
   import math
   import numpy as np

   num = math.sqrt(4*m*k - c**2);
   den = 2*m;
   y = 2*np.multiply(np.exp((-c/den)*t),np.cos((num/den)*t));
   J = np.inner((y-obs),((y-obs).T));
   return J




