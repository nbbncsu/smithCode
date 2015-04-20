%
%                            Example7_15.m
%
%
% Construct synthetic data and sampling distribution for the damping
% parameter C as illustrated in Example 7.15.
%
% Required function: spring_function.m
%

  clear all
  close all
  
%
% Specify true parameter and measurement variance values.
%
  m = 1;
  k = 20.5;
  c = 1.5;
  sigma = .1;
  var = sigma^2;

%
% Construct synthetic data, sensitivity value, and covariance value V.
%

  num = sqrt(4*m*k - c^2);
  den = 2*m;
  
  t = 0:.01:5;
  error = sigma*randn(size(t));
  n = length(t);

  y = 2*exp((-c/den)*t).*cos((num/den)*t);
  dydc = -(t/m).*exp((-c/den)*t).*cos((num/den)*t) + (c*t/(m*num)).*exp((-c/den)*t).*sin((num/den)*t);
  V = inv(dydc*dydc')*var;
  obs = y + error;
  
%
% Construct sampling distribution and determine optimal value of c
% for this data.

  C = [1.4:.001:1.6];
  sample_dist = normpdf(C,1.5,sqrt(V));

  modelfun = @(c) spring_function(c,m,k,obs,t);
  [c_opt,fval] = fminsearch(modelfun,c);

%
% Plot results.
%

  figure(1)
  plot(t,y,'k','linewidth',2)
  hold on
  plot(t,obs,'x','linewidth',1)
  plot(t,0*y,'k','linewidth',2)
  hold off
  set(gca,'Fontsize',[20]);
  legend('Model','Synthetic Data','Location','NorthEast')
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
  plot(C,sample_dist,'linewidth',2)
  set(gca,'Fontsize',[18]);
  xlabel('Damping Parameter c')
  ylabel('Sampling Density')

