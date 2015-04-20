% Jordan E. Massad
% Oct 18, 2004
% Phase Fraction Rate Law (all temperatures)
%
%RATEODE
%dx = rateODE(t,xfrac,mparams,cparams)   
%Input 
% t - time point/vector.  
% xfrac - phase fractions/temperature vector.  
% mparams - material parameters.  
% cparams - computed parameters.
%Output
% dx(1) - rate of M- production.
% dx(2) - rate of M+ production.
% dx(3) - rate of temperature change.
function dx = rateODE(t,xfrac,mparams,cparams)  

dx = zeros(3,1);     
Minv = zeros(3);
%Material Parameters
EM = mparams(2); Er = cparams(1);
eT = mparams(3); du = cparams(2);
h = mparams(5); pcv = mparams(6);
rhoA = mparams(9); rhoM = mparams(10); 
tau = cparams(10);  Tc = cparams(5);

T = xfrac(3); 
xM = xfrac(1)+xfrac(2);

%Phase-dependent Specific Heat
% pcv is for Martensite by default and cA-cM=dc.
TR = cparams(2)/cparams(3);
% dc = 0.0*mparams(6); 
dc = 0.29*mparams(6);  
pcv = pcv+dc - dc*xM;

%Scale input tensile stress to shear stress according to crystal orientation.
%Shift applied stress to effective stress.
sig = cparams(11)*input_sigma(t) - cparams(12);
% sig = cparams(11)*040 - cparams(12);

%External Temperature
Te = input_Te(t); 

%Radiative Transfer
omega = 2e-3*(1e3/1e1+1/1.778);         %Geometry Factor
emiss = 0.10*5.67051e-8;                %Emissivity
h = h + emiss*omega*(T^2+Te^2)*(T+Te);  %HTC [MW/m^3/K]

%Joule Heating
rhot = 1e-8*(rhoA + (rhoM-rhoA)*xM);
Joule = input_J(t)/rhot*1e-6;       %Voltage Input [MW/m^3]
% Joule = rhot*input_J(t);            %Current Input [MW/m^3]

%Total Heat Transfer
H = -h*(T-Te) + Joule;

%Transformation Likelihoods
[Pm,Pp,pAm,pAp] = pFuns(sig,T,mparams,cparams); 

%Transformation Enthalpies
Lm = (1-Er)/EM/2*sig^2 - eT*sig + du + dc*(T-TR); 
Lp= Lm + 2*eT*sig;

%ODE "Mass" Matrix Inverse
Minv(1,1) = 1/tau;
Minv(2,2) = Minv(1,1);
Minv(3,1) = Lm/(tau*pcv);
Minv(3,2) = Lp/(tau*pcv);
Minv(3,3) = 1/pcv;

% xA=1-xM;
dx(1,:) = (-(pAm+Pm)*xfrac(1) - pAm*xfrac(2) + pAm) + (T<=Tc).*(Pp*xfrac(2)-0).*(xM==1);
dx(2,:) = (-(pAp+Pp)*xfrac(2) - pAp*xfrac(1) + pAp) + (T<=Tc).*(Pm*xfrac(1)-0).*(xM==1);
dx(3,:) = H;

dx = Minv*dx;

% End rateODE.m