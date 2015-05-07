% Jordan E. Massad
% Oct 18, 2004
% Compute Single-crystal Quantities for Polycrystal Integrations 
%
%
%POLYFUN
%[ep,T,xm,xp] = polyFun(tspan,Tini,mparams,cparams)   
%Input 
% tspan - input time vector
% Tini - initial temperature point
% mparams - material parameters:  
% cparams - computed parameters.
% alpha - crystal orientation angle.
%Output
% ep - strain output vector
% T - internal temperature vector
% xm - M- phase fraction vector
% xp - M+ phase fraction vector
function [ep,T,xm,xp] = polyFun(tspan,Tini,mparams,cparams)   

% Solve for phase fractions:
dt = abs(tspan(end)-tspan(1));
opts =[];
% opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
opts = odeset(opts,'MaxStep',2e-2*dt);
% opts = odeset(opts,'InitialStep',2e-2*dt);
% opts = odeset('Stats','on');

%%Estimate Initial Conditions
[xplus0, xminus0] = phase_ini(tspan,Tini,mparams,cparams);

%%Solve ODE
warning off MATLAB:singularMatrix
[t,xfrac] = ode15s(@rateODE,tspan,[xminus0 xplus0 Tini],opts,mparams,cparams);  

%Scale tensile-to-shear and shift to effective stress:
sig = cparams(11)*input_sigma(t) + cparams(12);

%Compute Average Phase Strains
[epA,epMp,epMm] = mepsv(sig,xfrac(:,3),mparams,cparams);
%Compute Strain 
ep = epA + (epMm-epA).*xfrac(:,1) + (epMp-epA).*xfrac(:,2);
%Compute/Add Thermal Strains
lamA = mparams(7); lamM = mparams(8);
lamt = lamA + (lamM-lamA)*sum(xfrac(:,1:2)');
ep = ep + lamt(:).*(xfrac(:,3)-xfrac(1,3))*1e-6;

ep = ep(:)';
T = xfrac(:,3); T = T(:)';
xm = xfrac(:,1); xm = xm(:)';
xp = xfrac(:,2); xp = xp(:)';

% End polyFun.m