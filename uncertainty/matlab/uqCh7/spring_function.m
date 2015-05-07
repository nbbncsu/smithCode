  

  function J = spring_function(c,m,k,obs,t)

  num = sqrt(4*m*k - c^2);
  den = 2*m;
  y = 2*exp((-c/den)*t).*cos((num/den)*t);
  J = (y-obs)*((y-obs)');
