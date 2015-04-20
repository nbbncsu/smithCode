%
%                           Example12_4.m
%
%
% Code computes and plots the optimal parameters for the heat model (3.21)
% simultaneously using copper and aluminum data using the new boundary 
% condition (12.2).
%
% Required functions: heat_fun_newbc.m
% Required data: final_cu_data.txt
%                final_al_data.txt

  clear all
  close all

  global data_cu data_al xdata
  
%
% Load the data and construct the x datapoints.
%
  
  load final_cu_data.txt
  load final_al_data.txt

  data_cu = final_cu_data(2:16);
  data_al = final_al_data(2:16);
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  u_amb_cu = final_cu_data(17);
  u_amb_al = final_al_data(17);
  xvals = [0:.1:70];

%
% Input dimensions and material constants
%

  a = 0.95;      % cm
  b = 0.95;      % cm
  L = 70.0;      % cm
  k_cu = 4.01;   % W/cm C
  k_al = 2.37;   % W/cm C
  n = 15;        % Number of measurements
  p = 3;         % Number of parameters
  
  h_init = 0.0016;
  Tsource_init = -40.0;
  eta_init = -.0034;
  q_init = [h_init Tsource_init eta_init];

%
% Optimize parameters
%

  options = optimoptions('fmincon','Algorithm','active-set');
  modelfun = @(q)heat_fun_newbc(q,a,b,L,k_al,k_cu,u_amb_cu,u_amb_al);
  [q_opt,fval] = fminsearch(modelfun,q_init);


  h = q_opt(1)
  Tsource = q_opt(2)
  eta = q_opt(3)
  
%
% Compute solution using optimal parameter values.
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

  uvals_cu = c1_cu*exp(-gamma_cu*xvals) + c2_cu*exp(gamma_cu*xvals) + u_amb_cu;
  uvals_al = c1_al*exp(-gamma_al*xvals) + c2_al*exp(gamma_al*xvals) + u_amb_al;

  uvals_cu_data = c1_cu*exp(-gamma_cu*xdata) + c2_cu*exp(gamma_cu*xdata) + u_amb_cu;
  uvals_al_data = c1_al*exp(-gamma_al*xdata) + c2_al*exp(gamma_al*xdata) + u_amb_al;

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
  plot(xvals,uvals_al,'linewidth',2)
  axis([0 70 20 uvals_al(1)])
  hold on
  plot(xdata,data_al,'o','linewidth',5)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('Distance (cm)')
  ylabel('Temperature (^oC)')
  legend(' Model',' Data','Location','Northeast')
  
