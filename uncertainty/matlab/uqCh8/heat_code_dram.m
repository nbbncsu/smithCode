%
%                       heat_code_dram.m
%
%
% This code illustrates the implementation of the DRAM algorithm used 
% used to construct the chains and marginal densities in Figures 8.11
% and 8.12 for the heat Example 8.12.

%
% Required functions: kde.m
%                     heatss.m
%
% Required package: The DRAM package is available at the websites
%           https://wiki.helsinki.fi/display/inverse/Adaptive+MCMC
%           http://helios.fmi.fi/~lainema/mcmc/
%
% Required data: final_al_data.txt
%

  clear all
  close all

  load final_al_data.txt

  udata = final_al_data(2:16);
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  u_amb = final_al_data(17);

%
% Input dimensions and material constants
%

  a = 0.95;   % cm
  b = 0.95;   % cm
  L = 70.0;   % cm
  k = 2.37;   % W/cm C
  h = 0.00191;
  Q = -18.41; 
  n = 15;
  p = 2;

%
% Construct constants and solution
%

  gamma = sqrt(2*(a+b)*h/(a*b*k));
  gamma_h = (1/(2*h))*gamma;
  f1 = exp(gamma*L)*(h + k*gamma);
  f2 = exp(-gamma*L)*(h - k*gamma);
  f3 = f1/(f2 + f1);
  f1_h = exp(gamma*L)*(gamma_h*L*(h+k*gamma) + 1 + k*gamma_h);
  f2_h = exp(-gamma*L)*(-gamma_h*L*(h-k*gamma) + 1 - k*gamma_h);
  c1 = -Q*f3/(k*gamma);
  c2 = Q/(k*gamma) + c1;
  f4 = Q/(k*gamma*gamma);
  den2 = (f1+f2)^2;
  f3_h = (f1_h*(f1+f2) - f1*(f1_h+f2_h))/den2;
  c1_h = f4*gamma_h*f3 - (Q/(k*gamma))*f3_h;
  c2_h = -f4*gamma_h + c1_h;
  c1_Q = -(1/(k*gamma))*f3;
  c2_Q = (1/(k*gamma)) + c1_Q;

  uvals_data = c1*exp(-gamma*xdata) + c2*exp(gamma*xdata) + u_amb;
  uvals_Q_data = c1_Q*exp(-gamma*xdata) + c2_Q*exp(gamma*xdata);
  uvals_h_data = c1_h*exp(-gamma*xdata) + c2_h*exp(gamma*xdata) + gamma_h*xdata.*(-c1*exp(-gamma*xdata) + c2*exp(gamma*xdata));

  res = udata - uvals_data;

  sens_mat = [uvals_Q_data; uvals_h_data];
  sigma2 = (1/(n-p))*res*res';
  V = sigma2*inv(sens_mat*sens_mat');

%
% Construct the xdata, udata and covariance matrix V employed in DRAM.
% Set the options employed in DRAM.
%

  clear data model options

  data.xdata = xdata';
  data.ydata = udata';
  tcov = V;
  tmin = [Q; h];

  params = {
    {'q1',tmin(1), -20}
    {'q2',tmin(2), 0}};
  model.ssfun = @heatss;
  model.sigma2 = sigma2;
  options.qcov = tcov;
  options.nsimu = 10000;
  options.updatesigma = 1;
  N = 10000;

%
% Run DRAM to construct the chains (stored in chain) and measurement
% variance (stored in s2chain).
%

  [results,chain,s2chain] = mcmcrun(model,data,params,options);
  
%
% Construct the densities for Q and h.
%

  Qvals = chain(:,1);
  hvals = chain(:,2);
  range_Q = max(Qvals) - min(Qvals);
  range_h = max(hvals) - min(hvals);
  Q_min = min(Qvals)-range_Q/10;
  Q_max = max(Qvals)+range_Q/10;
  h_min = min(hvals)-range_h/10;
  h_max = max(hvals)+range_h/10;
  [bandwidth_Q,density_Q,Qmesh,cdf_Q]=kde(Qvals);
  [bandwidth_h,density_h,hmesh,cdf_h]=kde(hvals);

%
% Plot the results.
%

  figure(1); clf
  mcmcplot(chain,[],results,'chainpanel');

  figure(2); clf
  mcmcplot(chain,[],results,'pairs');

  cov(chain)
  chainstats(chain,results)

  figure(3)
  plot(Qvals,'-')
  set(gca,'Fontsize',[22]);
  axis([0 N -19.3 -17.5])
  xlabel('Chain Iteration')
  ylabel('Parameter Q')

  figure(4)
  plot(hvals,'-')
  set(gca,'Fontsize',[22]);
  axis([0 N 1.84e-3 2e-3])
  xlabel('Chain Iteration')
  ylabel('Parameter h')

  figure(5)
  plot(Qmesh,density_Q,'k-','linewidth',3)
  axis([-19.5 -17.5 0 3])
  set(gca,'Fontsize',[22]);
  xlabel('Parameter Q')

  figure(6)
  plot(hmesh,density_h,'k-','linewidth',3)
  axis([1.8e-3 2e-3 0 3e4])
  set(gca,'Fontsize',[22]);
  xlabel('Parameter h')

  figure(7)
  scatter(Qvals,hvals)
  box on
  axis([-19.2 -17.6 1.83e-3 2e-3])
  set(gca,'Fontsize',[23]);
  xlabel('Parameter Q')
  ylabel('Parameter h')




 


