
  function J = heat_fun_al(q,a,b,L,k_al,u_amb_al)

  global data_al xdata

  h = q(1);
  Q = q(2);

%
%  Construct constants and solution
%

  gamma_al = sqrt(2*(a+b)*h/(a*b*k_al));
  f1_al = exp(gamma_al*L)*(h + k_al*gamma_al);
  f2_al = exp(-gamma_al*L)*(h - k_al*gamma_al);
  f3_al = f1_al/(f2_al + f1_al);
  c1_al = -Q*f3_al/(k_al*gamma_al);
  c2_al = Q/(k_al*gamma_al) + c1_al;

  uvals_al = c1_al*exp(-gamma_al*xdata) + c2_al*exp(gamma_al*xdata) + u_amb_al;
  res_al = data_al - uvals_al;
  J = res_al*res_al';

