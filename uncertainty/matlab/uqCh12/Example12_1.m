%
%                           Example12_1.m
%
%
% Code computes and plots the optimal parameters for the heat model (3.21)
% using copper data.
%
% Required functions: heat_fun_cu.m
% Required data: final_cu_data.txt
%

  clear all
  close all

  global data_cu xdata
  
%
% Load the data and construct the x datapoints.
%
  
  load final_cu_data.txt

  data_cu = final_cu_data(2:16);
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  u_amb_cu = final_cu_data(17);
  xvals = [0:.1:70];

%
% Input dimensions and material constants
%

  a = 0.95;      % cm
  b = 0.95;      % cm
  L = 70.0;      % cm
  k_cu = 4.01;   % W/cm C
  n = 15;        % Number of measurements
  p = 2;         % Number of parameters
  
  h_init = 0.00183;
  Q_init = -15.93;
  q_init = [h_init Q_init];

%
% Optimize parameters
%

  modelfun = @(q)heat_fun_cu(q,a,b,L,k_cu,u_amb_cu);
  [q_opt,fval] = fminsearch(modelfun,q_init);

  q_opt 

  h = q_opt(1);
  Q = q_opt(2);

%
% Compute solution using optimal parameter values.
%

  gamma_cu = sqrt(2*(a+b)*h/(a*b*k_cu));
  f1_cu = exp(gamma_cu*L)*(h + k_cu*gamma_cu);
  f2_cu = exp(-gamma_cu*L)*(h - k_cu*gamma_cu);
  f3_cu = f1_cu/(f2_cu + f1_cu);
  c1_cu = -Q*f3_cu/(k_cu*gamma_cu);
  c2_cu = Q/(k_cu*gamma_cu) + c1_cu;

  uvals_cu = c1_cu*exp(-gamma_cu*xvals) + c2_cu*exp(gamma_cu*xvals) + u_amb_cu;
  uvals_cu_data = c1_cu*exp(-gamma_cu*xdata) + c2_cu*exp(gamma_cu*xdata) + u_amb_cu;

  res_cu = data_cu - uvals_cu_data;

%
% Plot the results
%

  figure(1)
  plot(xvals,uvals_cu,'linewidth',2)
  axis([0 70 20 uvals_cu(1)])
  hold on
  plot(xdata,data_cu,'o','linewidth',5)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('Distance (cm)')
  ylabel('Temperature (^oC)')
  legend(' Model',' Data','Location','Northeast')
  
  figure(2)
  plot(xdata,res_cu,'o','linewidth',6)
  axis([0 70 -.4 .45])
  hold on
  plot(xvals,0*xvals,'linewidth',2)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('Distance (cm)')
  ylabel('Residuals (^oC)')
  


