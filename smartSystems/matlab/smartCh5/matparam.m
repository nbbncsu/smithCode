% Jordan E. Massad
% Oct 18, 2004
% Define Model Material Parameters
%
%
%MATPARAM
%[mparams, cparams] = matparam(cs)   
%Input 
% n - null input.
%Output
% mparams - material parameters (see below)
function mparams = matparam(n)    

mparams = zeros(19,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material Parameters
mparams(1) = 37.2e3;                %EA [MPa]
mparams(2) = 25.5e3;                %EM [MPa]
mparams(3) = 0.0300;                %eT []
mparams(4) = 100.00;                %dsig (Mean Gap Stress) [MPa]
h = 1e0;                            %Convective Heat Transfer Coeff [W/m^2/K]
mparams(5) = 2e-3*(1e3/5e2+1/7)*h;  %Geometry-dependent HTC [MW/m^3/K]
                                	%   Formula: surf_area/volume * h.
rho = 6.45e3;                   	%SMA density [kg/m^3]
cv = 2.30e2;	                    %Specific heat (constant volume) [J/(kg*K)]
mparams(6) = rho*cv*1e-6;           %pcv [MJ/m^3/K]
mparams(7) = 11.0;                  %lamA [1e-6/K]  Austenite Thermal Expansion
mparams(8) = 6.6;               	%lamM [1e-6/K]  Martensite Thermal Expansion
mparams(9) = 100.0;             	%rhoA [1e-8 Ohm*m]  Austenite Resistivity
mparams(10) = 80.0;             	%rhoM [1e-8 Ohm*m]  Martensite Resistivity
mparams(11) = 248.56;           	%T_M [K]  Measured Mf<T_M<Ms
mparams(12) = 260;              	%T_A or Teq [K]  Measured As<T_A<Af or Ms<Teq<Af
mparams(13) = 274.0;            	%s_hi [MPa] Loading Transf. Point at T_hi>T_lo
mparams(14) = 293.0;             	%T_hi [K]  Temperature for s_hi
mparams(15) = 200.58;           	%s_lo [MPa] (scale measured tensile data to shear)
mparams(16) = 243.0;             	%T_lo [K]  Temperature for s_lo
mparams(17) = 52^2; 	            %Gap Variance [MPa^2] {52 default)
mparams(18) = 55^2;  	            %Effective Stress Variance [MPa^2] {45 default}

mparams(19) = 0.1*mparams(6);        %Specific Heat Difference (dc) [J/(kg*K)]


%  End of matparam.m