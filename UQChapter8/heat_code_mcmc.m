%
%                       heat_code_mcmc.m
%
%
% This code illustrates the implementation of the Metropolis algorithm 8.5
% used to construct the chains in Figure 8.10 for the heat Example 8.12.
% The code takes approximately 30 seconds to run.  This can be reduced
% by reducing the number of chain iterations.
%
% Required functions: kde.m
% Required data: final_al_data.txt
%

  clear all
  close all

  load final_al_data.txt

  data = final_al_data(2:16);
  xdata = [10 14 18 22 26 30 34 38 42 46 50 54 58 62 66];
  xvals = [10:.1:70];
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

  uvals = c1*exp(-gamma*xvals) + c2*exp(gamma*xvals) + u_amb;
  uvals_data = c1*exp(-gamma*xdata) + c2*exp(gamma*xdata) + u_amb;
  uvals_Q_data = c1_Q*exp(-gamma*xdata) + c2_Q*exp(gamma*xdata);
  uvals_h_data = c1_h*exp(-gamma*xdata) + c2_h*exp(gamma*xdata) + gamma_h*xdata.*(-c1*exp(-gamma*xdata) + c2*exp(gamma*xdata));

  res = data - uvals_data;

  sens_mat = [uvals_Q_data; uvals_h_data];
  sigma2 = (1/(n-p))*res*res';
  V = sigma2*inv(sens_mat*sens_mat');

%
% Set MCMC parameters
%
  N = 1e+5;
  R = chol(V);
  q_old = [Q;h];
  SS_old = res*res';
  n0 = 0.001;
  sigma02 = sigma2;
  aval = 0.5*(n0 + 15);
  bval = 0.5*(n0*sigma02 + SS_old);
  sigma2 = 1/gamrnd(aval,1/bval);
  accept = 0;

%
%  Run the Metropolis algorithm for N iterations.
%
  
  for i = 1:N
    z = randn(2,1); 
    q_new = q_old + R*z;
    Q = q_new(1,1);
    h = q_new(2,1);
    gamma = sqrt(2*(a+b)*h/(a*b*k));
    f1 = exp(gamma*L)*(h + k*gamma);
    f2 = exp(-gamma*L)*(h - k*gamma);
    f3 = f1/(f2 + f1);
    c1 = -Q*f3/(k*gamma);
    c2 = Q/(k*gamma) + c1;
    uvals_data = c1*exp(-gamma*xdata) + c2*exp(gamma*xdata) + u_amb;
    res = data - uvals_data;
    SS_new = res*res';
    u_alpha = rand(1);
    term = exp(-.5*(SS_new-SS_old)/sigma2);
    alpha = min(1,term);
    if u_alpha < alpha
      Q_MCMC(:,i) = [Q; h];
      q_old = q_new;
      SS_old = SS_new;
      accept = accept + 1;
    else
      Q_MCMC(:,i) = q_old;
    end
    Sigma2(i) = sigma2;
    bval = 0.5*(n0*sigma02 + SS_old);
    sigma2 = 1/gamrnd(aval,1/bval);
  end

  Qvals = Q_MCMC(1,:);
  hvals = Q_MCMC(2,:);
  
%
% Use kde to construct densities for Q and h.
%

  range_Q = max(Qvals) - min(Qvals);
  range_h = max(hvals) - min(hvals);
  Q_min = min(Qvals)-range_Q/10;
  Q_max = max(Qvals)+range_Q/10;
  h_min = min(hvals)-range_h/10;
  h_max = max(hvals)+range_h/10;
  [bandwidth_Q,density_Q,Qmesh,cdf_Q]=kde(Qvals);
  [bandwidth_h,density_h,hmesh,cdf_h]=kde(hvals);


  accept/N

%
% Plot solutions
%
 
  figure(1)
  plot(Qvals,'-','linewidth',2)
  set(gca,'Fontsize',[22]);
  axis([0 N -19.3 -17.5])
  xlabel('Chain Iteration')
  ylabel('Parameter Q')

  figure(2)
  plot(hvals,'-','linewidth',2)
  set(gca,'Fontsize',[22]);
  axis([0 N 1.84e-3 2e-3])
  xlabel('Chain Iteration')
  ylabel('Parameter h')

  figure(3)
  plot(Sigma2)
  set(gca,'Fontsize',[22]);
  title('Measurement Error Variance \sigma^2')

  figure(4)
  plot(Qmesh,density_Q,'k-','linewidth',3)
  axis([-19.5 -17.5 0 3])
  set(gca,'Fontsize',[22]);
  xlabel('Parameter Q')

  figure(5)
  plot(hmesh,density_h,'k-','linewidth',3)
  axis([1.8e-3 2e-3 0 3e4])
  set(gca,'Fontsize',[22]);
  xlabel('Parameter h')

  figure(6)
  scatter(Qvals,hvals)
  box on
  axis([-19.2 -17.6 1.83e-3 2e-3])
  set(gca,'Fontsize',[22]);
  xlabel('Parameter Q')
  ylabel('Parameter h')

