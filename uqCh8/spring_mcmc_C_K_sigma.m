%
%                       spring_mcmc_C_K_sigma.m
%
%
% This code illustrates the implementation of the Metropolis algorithm
% in Example 8.7, Case ii.  Here we are sampling the stiffness and damping
% parameter q = [C,K] as well as the variance sigma^2 for the measurement
% error. 
%

  clear all

%
% Input true stiffness and damping parameters and error variance
%

  K = 20.5;
  C = 1.5;
  sigma = .1;
  var = sigma^2;
  N = 10000;

  factor = sqrt(K - C^2/4);
  denom = sqrt(4*K - C^2);

%
% Construct solution on [0,5], the analytic sensitivity relations, and
% the observations, which include the error.
%

  t = 0:.01:5;
  n = length(t);

  y = 2*exp(-C*t/2).*cos(factor*t);
  dydc = exp(-C*t/2).*((C*t/denom).*sin(factor*t) - t.*cos(factor*t));
  dydk = (-2*t/denom).*exp(-C*t/2).*sin(factor*t);

  error = sigma*randn(size(t));
  obs = y + error;

%
% Construct the sensitivity matrix chi and covariance matrix V.
%

  chi = [dydc' dydk'];
  V = inv(chi'*chi)*var;
  R = chol(V);

  q_old = [C;K];
  SS_old = (obs-y)*(obs-y)';

%
% Construct the parameters used when sampling the error variance in the Metropolis algorithm.
%

  n0 = .001;
  sigma02 = var;
  aval = 0.5*(n0 + n);
  bval = 0.5*(n0*sigma02 + SS_old);
  sigma2 = 1/gamrnd(aval,1/bval);

%
% Construct a Metropolis chain of length N
%

  for i = 1:N
    z = randn(2,1); 
    q_new = q_old + R*z;
    C_new = q_new(1,1);
    K_new = q_new(2,1);
    factor = sqrt(K_new - C_new^2/4);
    y_new = 2*exp(-C_new*t/2).*cos(factor*t);
    u_alpha = rand(1);
    SS_new = (obs-y_new)*(obs-y_new)';
    term = exp(-.5*(SS_new-SS_old)/sigma2);
    alpha = min(1,term);
    if u_alpha < alpha
      Q_MCMC(:,i) = [C_new; K_new];
      q_old = q_new;
      SS_old = SS_new;
    else
      Q_MCMC(:,i) = q_old;
    end
    Sigma2(i) = sigma2;
    bval = 0.5*(n0*sigma02 + SS_old);
    sigma2 = 1/gamrnd(aval,1/bval);
  end

  Cvals = Q_MCMC(1,:);
  Kvals = Q_MCMC(2,:);
  
%
% Use kde to construct densities for C and K
%

  range_C = max(Cvals) - min(Cvals);
  range_K = max(Kvals) - min(Kvals);
  C_min = min(Cvals)-range_C/10;
  C_max = max(Cvals)+range_C/10;
  K_min = min(Kvals)-range_K/10;
  K_max = max(Kvals)+range_K/10;
  [bandwidth_C,density_C,Cmesh,cdf_C]=kde(Cvals);
  [bandwidth_K,density_K,Kmesh,cdf_K]=kde(Kvals);

%
% Plot results.  Note that the axis limits may change as a function of
% each new generated data set.
%

  figure(1)
  plot(Cvals,'-','linewidth',1)
  set(gca,'Fontsize',[22]);
  axis([0 N 1.43 1.57])
  xlabel('Chain Iteration')
  ylabel('Damping Parameter C')

  figure(2)
  plot(Kvals,'-','linewidth',1)
  set(gca,'Fontsize',[22]);
%  axis([0 N 20.2 20.9])
  xlabel('Chain Iteration')
  ylabel('Stiffness Parameter K')

  figure(3)
  plot(Sigma2)
  set(gca,'Fontsize',[22]);
  title('Measurement Error Variance \sigma^2')

  figure(4)
  plot(Cmesh,density_C,'k-','linewidth',3)
  axis([1.4 1.6 0 25])
  set(gca,'Fontsize',[22]);
  xlabel('Damping Parameter C')

  figure(5)
  plot(Kmesh,density_K,'k-','linewidth',3)
%  axis([20 21 0 5])
  set(gca,'Fontsize',[22]);
  xlabel('Stiffness Parameter K')

  figure(6)
  scatter(Cvals,Kvals)
  box on
 % axis([-19.2 -17.6 1.83e-3 2e-3])
  set(gca,'Fontsize',[22]);
  xlabel('Damping Parameter C')
  ylabel('Stiffness Parameter K')

