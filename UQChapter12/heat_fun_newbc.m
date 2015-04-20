
  function J = heat_fun_newbc(q,a,b,L,k_al,k_cu,u_amb_cu,u_amb_al)

  global data_cu data_al xdata

  h = q(1);
  Tsource = q(2);
  eta = q(3);

%
%  Construct constants and solution
%

  gamma_cu = sqrt(2*(a+b)*h/(a*b*k_cu));
  Gamma_cu = (k_cu*gamma_cu + eta)/(k_cu*gamma_cu - eta);
  f1_cu = exp(gamma_cu*L)*Gamma_cu*(h + k_cu*gamma_cu);
  f2_cu = exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
  f3_cu = exp(gamma_cu*L)*(h + k_cu*gamma_cu);
  f4_cu = f3_cu/(f2_cu + f1_cu);
  c1_cu = f4_cu*eta*(u_amb_cu - Tsource)/(eta - k_cu*gamma_cu);
  c2_cu = eta*(Tsource - u_amb_cu)/(eta - k_cu*gamma_cu) + c1_cu*Gamma_cu;

  gamma_al = sqrt(2*(a+b)*h/(a*b*k_al));
  Gamma_al = (k_al*gamma_al + eta)/(k_al*gamma_al - eta);
  f1_al = exp(gamma_al*L)*Gamma_al*(h + k_al*gamma_al);
  f2_al = exp(-gamma_al*L)*(h - k_al*gamma_al);
  f3_al = exp(gamma_al*L)*(h + k_al*gamma_al);
  f4_al = f3_al/(f2_al + f1_al);
  c1_al = f4_al*eta*(u_amb_al - Tsource)/(eta - k_al*gamma_al);
  c2_al = eta*(Tsource - u_amb_al)/(eta - k_al*gamma_al) + c1_al*Gamma_al;

  uvals_cu = c1_cu*exp(-gamma_cu*xdata) + c2_cu*exp(gamma_cu*xdata) + u_amb_cu;
  uvals_al = c1_al*exp(-gamma_al*xdata) + c2_al*exp(gamma_al*xdata) + u_amb_al;

  res_cu = data_cu - uvals_cu;
  res_al = data_al - uvals_al;
  J = res_cu*res_cu' + res_al*res_al';

