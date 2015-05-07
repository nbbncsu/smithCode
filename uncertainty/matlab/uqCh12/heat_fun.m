
  function J = heat_fun(q,a,b,L,k_al,k_cu,u_amb_al,u_amb_cu)

  global data_cu data_al xdata

  h = q(1);
  Q = q(2);

%
%  Construct constants and solution
%

  gamma_cu = sqrt(2*(a+b)*h/(a*b*k_cu));
  f1_cu = exp(gamma_cu*L)*(h + k_cu*gamma_cu);
  f2_cu = exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
  f3_cu = f1_cu/(f2_cu + f1_cu);
  c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
  c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

  gamma_al = sqrt(2*(a+b)*h/(a*b*k_al));
  f1_al = exp(gamma_al*L)*(h + k_al*gamma_al);
  f2_al = exp(-gamma_al*L)*(h - k_al*gamma_al);
  f3_al = f1_al/(f2_al + f1_al);
  c1_al = -Q*f3_al/(k_al*gamma_al);
  c2_al = Q/(k_al*gamma_al) + c1_al; 

  uvals_cu = c1_cu*exp(-gamma_cu*xdata) + c2_cu*exp(gamma_cu*xdata) + u_amb_cu;
  uvals_al = c1_al*exp(-gamma_al*xdata) + c2_al*exp(gamma_al*xdata) + u_amb_al;

  res_cu = data_cu - uvals_cu;
  res_al = data_al - uvals_al;
  J = res_cu*res_cu' + res_al*res_al';
