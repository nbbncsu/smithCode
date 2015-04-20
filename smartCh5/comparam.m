% Jordan E. Massad
% Oct 18, 2004
% Computed Model Parameters
%
%
%COMPARAM
%[cparams] = comparam(mparams)   
%Input 
% mparams - material parameters (see matparam.m).
%Output
% cparams - computed parameters (see below).
function [cparams] = comparam(mparams)    

cparams = zeros(12,1);

eT=mparams(3); dsig=mparams(4); 
T_M=mparams(11); T_A=mparams(12);
s_lo=mparams(15); s_hi=mparams(13);
T_hi=mparams(14); T_lo=mparams(16); 
Teq = mparams(12); TR = mparams(12);
dc = mparams(19); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computed Parameters
cparams(1) = mparams(2)/mparams(1);                             %Er    []
Er = cparams(1);

%%Identification Algorithms
%%%
%%%2 Loading Points Formulation
%%%
% cparams(2) = (Er-1)*(T_hi*s_lo^2-T_lo*s_hi^2) + ...
%     ((1-Er)*dsig-2*mparams(2)*eT)*(T_hi*s_lo-T_lo*s_hi);
% cparams(2) = cparams(2)/mparams(2)/(T_hi-T_lo)/2 + dsig*eT/2;  %du [MJ/m^3] Reference Internal Energy Difference
% 
% cparams(3) = (s_hi-s_lo)*((1-Er)*(s_lo+s_hi-dsig) + 2*mparams(2)*eT);
% cparams(3) = cparams(3)/mparams(2)/(T_hi-T_lo)/2;              %ds [MJ/m^3/K] Reference Entropy Difference

%%%
%%%2 Stress-dependent T_M Formulation
%%%
% Flo = s_lo*(eT + (1-Er)/2*(s_lo-dsig)/mparams(2) );
% Fhi = s_hi*(eT + (1-Er)/2*(s_hi-dsig)/mparams(2) );
% cparams(2) = (Fhi*T_lo-Flo*T_hi)/(T_hi-T_lo) + dsig*eT/2;      %du [MJ/m^3] Reference Internal Energy Difference 
% cparams(3) = (Fhi-Flo)/(T_hi-T_lo);                            %ds [MJ/m^3/K] Reference Entropy Difference

%%%
%%%Mf<T_M<Ms + As<T_A<Af Formulation
%%%
% cparams(2) = dsig*eT*(T_M+T_A)/(T_A-T_M)/2;                    %du [MJ/m^3] Reference Internal Energy Difference
% cparams(3) = dsig*eT/(T_A-T_M);                                %ds [MJ/m^3/K] Reference Entropy Difference
% T_hi=T_A; Teq=T_M;

%%%
%%%Mf<T_M<Ms + 1 Loading Point Formulation
%%%
% cparams(3) = 2*mparams(2)*eT + (1-Er)*(s_hi-dsig);
% cparams(3) = s_hi*cparams(3)/mparams(2)/(T_hi-T_M)/2;    		%ds [MJ/m^3/K] Reference Entropy Difference
% cparams(2) = T_M*cparams(3) + dsig*eT/2;                 		%du [MJ/m^3] Reference Internal Energy Difference

%%%
%%%Teq=mean(Ms,Af) + 1 Loading Point (at T_hi>Teq) Formulation
%%%
% Calculate Energy Parameters (for dc=0) 
cparams(3) = s_hi*(2*mparams(2)*eT + (1-Er)*(s_hi-dsig));
cparams(3) = (cparams(3)/mparams(2) - dsig*eT)/(T_hi-Teq)/2;   	%ds [MJ/m^3/K] Reference Entropy Difference
cparams(2) = Teq*cparams(3);                 		            %du [MJ/m^3] Reference Internal Energy Difference


% Nonzero dc Corrections
cparams(2) = cparams(2) + dc*(Teq*T_hi*log(Teq/T_hi)/(T_hi-Teq) + TR);
cparams(3) = cparams(3) + dc*((Teq*log(Teq/TR)-T_hi*log(T_hi/TR))/(T_hi-Teq) + 1);

%%%Calculate Transition Temperatures (nonzero dc)
opts=optimset('fzero'); opts=optimset(opts,'disp','off');

%%Tmin [K] (wM=0 lowest allowable Helmholtz temp)
dBTmin = mparams(2)*eT^2/2;
Tmin = fzero(@dBfun,TR,opts,cparams(2),cparams(3),dc,TR,dBTmin);
cparams(4) = max(1e-2,Tmin);                                    
%Alternative Computation
% T0 = Teq*exp(-cparams(3)/dc)/2
% Tmin = fzero(@dBfun,T0,opts,cparams(2),cparams(3),dc,TR,0);
% cparams(4) = Tmin*(Tmin<Teq);                                   

%%T_M [K] (sigA=0 2-well to 3-well temperature)
Tm = fzero(@dBfun,TR,opts,cparams(2),cparams(3),dc,TR,dsig*eT/2);
cparams(5) = max(2e-2,Tm);                                 		

%%Teq [K] (dB=0 equilibrium temperature)
Teq = fzero(@dBfun,TR,opts,cparams(2),cparams(3),dc,TR,0);		

%%T_A [K] (wM=eT superelastic boundary)
T_A = fzero(@dBfun,TR,opts,cparams(2),cparams(3),dc,TR,-dsig*eT/2);
cparams(6) = T_A;               	                            	

%%T_Asup [K] (wA=eT high temp beyond superelasticity)
dB3=eT/2/cparams(1)*(dsig - ((mparams(1)+mparams(2))*eT));
T_Asup = fzero(@dBfun,TR,opts,cparams(2),cparams(3),dc,TR,dB3);
cparams(7) = T_Asup;               	                            	

%Used in Boltzman integrals
kB = 1.3806503e-23;                                         	%Boltzman constant [J/K]
V =  1e-0*kB;                                             		%Lattice layer volume [m^3]
mhat = 3e3*V^(1/3);                            	                %Mass of a layer [kg]

cparams(8) = sqrt((V/kB)/mparams(1)/2)*1e3;               		%qA [sqrt(K)/MPa]
cparams(9) = sqrt((V/kB)/mparams(2)/2)*1e3;               		%qM [sqrt(K)/MPa]
cparams(10) = pi*sqrt(mhat/mparams(1)/V^(1/3))*1e-3;           	%Relaxation tau [s]

cparams(11) = 1;                                        		%Tensile-to-shear scaling factor = sin(2*alpha)/2 []
cparams(12) = 0;                                        		%Effective stress distributional shift [Pa]
                                                        		%  (values set in integration codes)

% End comparam()

%%Chemical Free Energy Difference
function dB = dBfun(T,du,ds,dc,TR,E)
dB = du -T*ds +dc*(T-TR-T.*log(T/TR)) - E;

%  End of comparam.m