  function J = heat_fun_cu_quad(q,a,b,L,k_cu,u_amb_cu)

%
% Function employed in the optimization routine for Example 12.6.
%
  global data_cu xdata

  h = q(1);
  Q = q(2);
  b0 = q(3);
  b1 = q(4);
  b2 = q(5);

%
%  Construct constants and solution
%

  gamma_cu = sqrt(2*(a+b)*h/(a*b*k_cu));
  f1_cu = exp(gamma_cu*L)*(h + k_cu*gamma_cu);
  f2_cu = exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
  f3_cu = f1_cu/(f2_cu + f1_cu);
  c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
  c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

  uvals_cu = c1_cu*exp(-gamma_cu*xdata) + c2_cu*exp(gamma_cu*xdata) + u_amb_cu + b0 + b1*xdata + b2*(xdata.^2); 
  res_cu = data_cu - uvals_cu;
  J = res_cu*res_cu';

