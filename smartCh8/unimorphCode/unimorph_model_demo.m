%
%      Unimorph_model_demo.m
%

clear all
close all
format short e

% Input the physical parameters and discretization limits. 

load 25volt.dat
load 50volt.dat
load 75volt.dat
load 100volt.dat

N = 16;
[ell,h_a,rho,Y,cD,gamma,Kb,YI,mil] = parameters_unimorph(N);
h = ell/N;
fval = 1;

V25 = h_a*1e+5*X25volt(:,1);
Dis25 = mil*X25volt(:,2);
V50 = h_a*1e+5*X50volt(:,1);
Dis50 = mil*X50volt(:,2);
V75 = h_a*1e+5*X75volt(:,1);
Dis75 = mil*X75volt(:,2);
V100 = h_a*1e+5*X100volt(:,1);
Dis100 = mil*X100volt(:,2);

%  Construct the system matrices and vectors.

[M,K,C] = matrix_construct_unimorph(N,ell,rho,Y,cD,gamma,YI);
F = force_unimorph(N,ell,fval);
B = b_input_unimorph(N,ell,Kb);
 
A_mat = [zeros(size(K)) eye(size(K))
         - inv(M)*K     -inv(M)*C];

F_vec = [zeros(size(F))
         inv(M)*F];

B_vec = [zeros(size(F))
         inv(M)*B];

% Integrate the system using the stiff ODE solver ode15s.m.
% The approximate solution is calculated at the point x on the beam.

t0 = 0;                             % Starting Time
tf = 4;                             % Final Time
nsteps = 250;                       % Number of time steps in interval
hsteps = (tf-t0)/nsteps;
y0 = zeros(2*(N+1),1);
x = ell;                            % Observation point on beam
bvn = solution_unimorph(x,N,ell);   % bvn contains the cubic splines evaluated at x
tspan = t0:hsteps:tf;

options = odeset('JConstant','on','Jacobian','on','BDF','on','RelTol',1e-9, 'AbsTol',1e-9);

[t,y25] = ode15s('yprime_unimorph',tspan,y0,options,A_mat,F_vec,B_vec,V25);
[t,y50] = ode15s('yprime_unimorph',tspan,y0,options,A_mat,F_vec,B_vec,V50);
[t,y75] = ode15s('yprime_unimorph',tspan,y0,options,A_mat,F_vec,B_vec,V75);
[t,y100] = ode15s('yprime_unimorph',tspan,y0,options,A_mat,F_vec,B_vec,V100);

yt25 = y25(:,1:N+1)';
dis25 = bvn*yt25;
volt25 = -max(V25)*sin(2*pi*t);
yt50 = y50(:,1:N+1)';
dis50 = bvn*yt50;
volt50 = -max(V50)*sin(2*pi*t);
yt75 = y75(:,1:N+1)';
dis75 = bvn*yt75;
volt75 = -max(V75)*sin(2*pi*t);
yt100 = y100(:,1:N+1)';
dis100 = bvn*yt100;
volt100 = -max(V100)*sin(2*pi*t);

Vpeak = 2*[max(V25) max(V50) max(V75) max(V100)];
Dispeak = [max(Dis25)-min(Dis25) max(Dis50)-min(Dis50) max(Dis75)-min(Dis75) max(Dis100)-min(Dis100)];
dispeak = [max(dis25)-min(dis25) max(dis50)-min(dis50) max(dis75)-min(dis75) max(dis100)-min(dis100)];
[P,S] = polyfit(Vpeak,Dispeak,1);
  
x1 = 120;
y1 = 1.1e-4;

figure(1)
subplot(2,2,1)
plot(volt25(51:251),dis25(51:251))
axis([-x1 x1 -y1 y1])
hold on
plot(V25(1:250),Dis25(1:250),'--')
plot([-x1 x1],[0 0],'k-',[0 0],[-y1 y1],'-')
xlabel('Voltage (V)')
ylabel('Tip Displacement (m)')
legend('Model','Data',1)
title('25 V Inputs')
hold off

subplot(2,2,2)
plot(volt50(51:251),dis50(51:251))
axis([-x1 x1 -y1 y1])
hold on
plot(V50(1:250),Dis50(1:250),'--')
plot([-x1 x1],[0 0],'k-',[0 0],[-y1 y1],'-')
xlabel('Voltage (V)')
ylabel('Tip Displacement (m)')
legend('Model','Data',1)
title('50 V Inputs')
hold off

subplot(2,2,3)
plot(volt75(51:251),dis75(51:251))
axis([-x1 x1 -y1 y1])
hold on
plot(V75(1:250),Dis75(1:250),'--')
plot([-x1 x1],[0 0],'k-',[0 0],[-y1 y1],'-')
xlabel('Voltage (V)')
ylabel('Tip Displacement (m)')
legend('Model','Data',1)
title('75 V Inputs')
hold off

subplot(2,2,4)
plot(volt100(51:251),dis100(51:251))
axis([-x1 x1 -y1 y1])
hold on
plot(V100(1:250),Dis100(1:250),'--')
plot([-x1 x1],[0 0],'k-',[0 0],[-y1 y1],'-')
h = gca;
get(h);
xlabel('Voltage (V)')
ylabel('Tip Displacement (m)')
legend('Model','Data',1)
title('100 V Inputs')
hold off


figure(2)
plot(Vpeak,dispeak)
hold on
plot(Vpeak,Dispeak,'x')
xlabel('Voltage (V)')
ylabel('Peak to Peak Tip Displacement (m)')
legend('Model','Data',4)
hold off 

