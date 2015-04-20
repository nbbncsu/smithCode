%                       heatss

  function ss = heatss(params,data)

  udata = data.ydata;
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66]';

% Input dimensions and material constants

  a = 0.95;   % cm
  b = 0.95;   % cm
  L = 70.0;   % cm
  k = 2.37;   % W/cm C
  Q = params(1); 
  h = params(2);
  n = 15;
  p = 2;
  u_amb = 21.2897; 

% Construct constants and solution

  gamma = sqrt(2*(a+b)*h/(a*b*k));
  f1 = exp(gamma*L)*(h + k*gamma);
  f2 = exp(-gamma*L)*(h - k*gamma);
  f3 = f1/(f2 + f1);
  c1 = -Q*f3/(k*gamma);
  c2 = Q/(k*gamma) + c1;

  uvals_data = c1*exp(-gamma*xdata) + c2*exp(gamma*xdata) + u_amb;

  res = udata - uvals_data;
  ss = res'*res;

