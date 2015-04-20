
%
%      Plate_model_demo.m
%

clear all
close all
format short e

%  Input the physical parameters and discretization limits.

[ellx,elly,actuator,act_sense,aval,kval,Ms,lams,alpha,cval,rho_p,E_p,...
    cD_p,gamma,nu,h,Kb,Nx,Ny,Nqx,Nqy,N_act,N_newt,N,bc] = parameters_plate;

t0 = 0;                             
t1 = 2.5;                        
nsteps1 = round((t1-t0)*100);  
hsteps1 = (t1-t0)/nsteps1;

Nxx = 20;
Nyy = 20;
obs = 441;
gx = 0:ellx/Nxx:ellx;
gy = 0:elly/Nyy:elly;
obs = 441;

freq_set = [10 12 19 26 37];

% Construct the mass matrix M, stiffness matrix K and damping matrix C.
% Also construct the disturbance vector F and the control input vector B.

[M,K,C] = matrix_construct_plate(Nx,Nqx,ellx,Ny,Nqy,elly,...
          rho_p,E_p,cD_p,gamma,nu,h,bc);

F = 10*force_plate(Nx,Nqx,ellx,Ny,Nqy,elly,rho_p,E_p,cD_p,gamma,nu,h,bc);

Bvn = solution_plate(gx',Nx,ellx,gy',Ny,elly,bc);

A_mat = [zeros(size(K)) eye(size(K))
         -inv(M)*K     -inv(M)*C];

F_vec = [zeros(size(F))
         inv(M)*F];

%
%  Integrate system using the stiff ODE solver ode15s.m. 
%
 
y0 = zeros(2*N,1);
tspan = t0:hsteps1:t1;  
options = odeset('JConstant','on','Jacobian','on','BDF','on','RelTol',1e-4, 'AbsTol',1e-4);
[t,y] = ode15s('yprime_plate',tspan,y0,options,A_mat,F_vec,freq_set);
yt1 = y(:,1:N);
ysolution = Bvn*yt1';

figure(1)
plot(t,ysolution(441,:))
xlabel('Time (s)')
ylabel('Displacement (m)')




