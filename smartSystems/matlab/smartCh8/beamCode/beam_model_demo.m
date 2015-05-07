%
%                Beam_model_demo.m
%

clear all
close all
format short e

% Input the physical parameters and discretization limits.  For the cantilever
% beam a total of N+1 = 17 cubic splines are used in the discretization.  

[ell,s,pe1,pe2,rho_b,rho_p,EI_b,EI_p,cD_b,cD_p,gamma_b,Kb,N,Nq] = parameters_beam;

%  Construct the mass matrix M, stiffness matrix K and damping
%  matrix C.   Also construct the disturbance vector F and the
%  control input vector B.  The disturbance was assumed to be constant
%  over the length of the beam (e.g., pressure) and periodic
%  in time with a frequency of 5 Hz.  This is applied for T=.74 sec
%  and then is cut to allow the beam vibrations to begin decay.
%
%  The second-order mass, stiffness, damping and input components
%  are then used to construct system matrices for the first-order
%  control system 
%
%     y'(t) = A_mat y(t) + B_vec u(t) + F_vec(t)

[M,K,C] = matrix_construct_beam(N,Nq,s,pe1,pe2,ell,...
    rho_b,rho_p,EI_b,EI_p,cD_b,cD_p,gamma_b);
F = force_beam(N,Nq,ell);
B = control_input_beam(N,Nq,ell,pe1,pe2,Kb);
 
A_mat = [zeros(size(K)) eye(size(K))
         - inv(M)*K     -inv(M)*C];

F_vec = [zeros(size(F))
         inv(M)*F];

B_vec = [zeros(size(F))
         inv(M)*B];
  
%  Integrate the system using the stiff ODE solver ode15s.m.
%  The approximate solution is calculated at the point x on the beam.

t0 = 0;                             % Starting Time
tf = 2.5;                           % Final Time
nsteps = 250;                       % Number of time steps in interval
hsteps = (tf-t0)/nsteps;
y0 = zeros(2*(N+1),1);
x = 3*ell/5;                        % Observation point on beam
bvn = solution_beam(x,N,ell);       % bvn contains the cubic splines evaluated at x=3*ell/5
tspan = t0:hsteps:tf;

options = odeset('reltol',1e-3,'abstol',1e-6,'vectorized','on','Jacobian',A_mat);
[t,y] = ode15s(@yprime_beam,tspan,y0,options,A_mat,F_vec);

yt = y(:,1:N+1);
w_ol = yt*bvn';
   
figure(1)
plot(t,w_ol)
h = gca;
get(h);
set(h,'FontSize',[16]);
xlabel('Time (s)');
ylabel('Displacement (m)');
