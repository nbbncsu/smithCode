function [ell,s,pe1,pe2,rho_b,rho_p,EI_b,EI_p,...
        cD_b,cD_p,gamma_b,Kb,N,Nq] = parameters_beam

%  The following dimensions and parameters are consistent with those of
%  the experimental beam B1TAPIH reported in [Banks, Wang, Inman, Slater,
%  Control-Theory and Advanced Technology, 10(4), 1994, pp. 873-892]. Metric
%  units are used in all cases.

h = .0016;                     % (m)
h_pe = .000254;                % (m)
b = .0203;                     % (m)
ell = .4573;                   % (m)

rho_b = .093;                  % (kg/m)
rho_p = .433;
EI_b = .491;                   % (N m^2)
EI_p = .793;
cD_b = .649e-5;                % (s N m^2)
cD_p = 1.255e-5;
gamma_b = .013;                % (s N/m^2)
Kb = 1.746e-2;                 % (N m/V)
 
s = 1;                         % Number of Patches
pe1 = .15;                     % Left edge of patch
pe2 = .25;                     % Right edge of patch

N = 16;
Nq = 16;