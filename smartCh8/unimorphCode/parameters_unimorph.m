function [ell,h_a,rho,Y,cD,gamma,Kb,YI,mil] = parameters(N)

mil = 25.4e-6;                 % (m)
ell = 0.03;                    % (m)
h_a = 52e-6;                   % (m)
h_i = 137e-6;                  % (m)
h = h_a + h_i;                 % (m)
h_pe = .000254;                % (m)
b = 0.013;                     % (m)

d_31 = 20e-12;                 % (C/N)
rho_a = 1.78e+3;
rho_i = 1.30e+3;
Y_a = 2.0e+9;
Y_i = 2.7e+9;
cD = 2.2848e-7;                % (s N m^2)
gamma = .005;                  % (s N/m^2)
Kb = 1.746e-2;                 % (N m/V)

s = 1;                         % Number of Patches
pe1 = 0;                       % Left edge of patch
pe2 = ell;                     % Right edge of patch

rho = (rho_a*h_a + rho_i*h_i)*b;
Zns = (Y_a*h_a^2-Y_i*h_i^2)/(2*(Y_i*h_i+Y_a*h_a));
M_kb= (Y_i*(-Zns^3+(h_i+Zns)^3)+Y_a*((h_a-Zns)^3+Zns^3));
Y = M_kb/((h_a-Zns)^3+(h_i+Zns)^3);
Kb = (1/2)*b*Y_a*d_31*(h_a-2*Zns); 
D_A = h_a - Zns;
D_I = h_i + Zns;
YI = Y*b*(D_A^3 + D_I^3)/3;
