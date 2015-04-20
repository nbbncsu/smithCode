
%                       spring_mcmc_C.m
%
%
% This code illustrates the implementation of the Metropolis algorithm
% in Example 8.7, Case i.
%

  clear all
  close all

%
% Input true stiffness and damping parameters and error variance
%

  K = 20.5;
  C = 1.5;
  C_old = C;
  C0 = C;
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
  y0 = y;
  dydc = exp(-C*t/2).*((C*t/denom).*sin(factor*t) - t.*cos(factor*t));

  error = sigma*randn(size(t));
  obs = y + error;

%
% Construct the sensitivity matrix chi and covariance matrix V.
%

  chi = [dydc'];
  V = inv(chi'*chi)*var;
  R = chol(V);

  SS_old = (obs-y)*(obs-y)';

%
% Construct a Metropolis chain of lenth N
%

  for i = 1:N
    z = randn(1);
    C_new = C_old + R*z;
    factor = sqrt(K - C_new^2/4);
    y_new = 2*exp(-C_new*t/2).*cos(factor*t);
    u_alpha = rand(1);
    SS_new = (obs-y_new)*(obs-y_new)';
    term = exp(-.5*(SS_new-SS_old)/var);
    alpha = min(1,term);
    if u_alpha < alpha
      Cvals(i) = C_new;
      C_old = C_new;
      SS_old = SS_new;
    else
      Cvals(i) = C_old;
    end
  end

%
% Use kde to construct the density for C from the Metropolis chain
%

  range_C = max(Cvals) - min(Cvals);
  C_min = min(Cvals)-range_C/10;
  C_max = max(Cvals)+range_C/10;
  [bandwidth_C,density_C,Cmesh,cdf_C]=kde(Cvals);

%
% Construct the posterior density directly by approximating Bayes' rule using a
% a midpoint quadrature rule.
%

  h = .001;
  Nq = .2/h;
  for i = 1:Nq
    C = 1.4 + (i-1)*h;
    factor = sqrt(K - C^2/4); 
    y = 2*exp(-C*t/2).*cos(factor*t);
    SS = (obs-y)*(obs-y)';
    for j = 1:Nq
      Cq = 1.4 + (j-1)*h;
      factor = sqrt(K - Cq^2/4); 
      yq = 2*exp(-Cq*t/2).*cos(factor*t);
      SSq = (obs-yq)*(obs-yq)';
      quadpt(j) = exp(-.5*(SSq-SS)/var);
    end
    posterior(i) = 1/(h*sum(quadpt));
    C_post_axis(i) = C;
  end

%
% Construct the asymptotic normal sampling distribution.
%

  mean = h*sum(C_post_axis.*posterior);
  sample_dist = normpdf(C_post_axis,C0,sqrt(V));

%
% Plot results. 
%

  figure(1)
  plot(t,y0,t,0*y0,'k','linewidth',2)
  hold on
  plot(t,obs,'x','linewidth',1)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('Time (s)')
  ylabel('Displacement')

  figure(2)
  plot(t,error,'x',t,0*y,'k','linewidth',2)
  hold on
  plot(t,2*sigma*ones(size(t)),'--r',t,-2*sigma*ones(size(t)),'--r','linewidth',2)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('Time (s)')
  ylabel('Residuals')
 
  figure(3)
  plot(Cvals,'-','linewidth',1)
  set(gca,'Fontsize',[22]);
  axis([0 N 1.43 1.57])
  xlabel('Chain Iteration')
  ylabel('Damping Parameter C')


  figure(4)
  hold on
  plot(C_post_axis,posterior,'b-','linewidth',3)
  plot(Cmesh,density_C,'g-.','linewidth',3) 
  plot(C_post_axis,sample_dist,'r--','linewidth',3) 
  hold off 
  axis([1.4 1.6 0 25])
  set(gca,'Fontsize',[20]);
  legend('Bayes','MCMC','OLS','Location','Northeast')
  xlabel('Damping Parameter C')


