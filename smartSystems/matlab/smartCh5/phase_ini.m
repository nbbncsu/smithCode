% Jordan E. Massad
% Oct 18, 2004
% Initial Phase Fraction Computation
%
%
%PHASE_INI
%[xplus0, xminus0] = phase_ini(tspan,Tini,mparams,cparams)   
%Input 
% Tini - initial environment temperature [K].
% mparams - material parameters.  
% cparams - computed parameters.  
%Output
% xplus0,xminus0 - initial M+/M- phase fraction.
function [xplus0, xminus0] = phase_ini(tspan,Tini,mparams,cparams)   

xplus0=0.0; xminus0=0.0;
T_M = cparams(5);
T_A = cparams(6);

%Initial Input
sig0 = input_sigma(tspan(1));
Te = input_Te(tspan(1:2));

%Transition Stresses
sigA = sigAfun(Tini,mparams,cparams);
sigM = sigA-mparams(4);
%Stress-dependent Transition Temperatures
if(sig0~=0)
	[T_Msig, T_Asig] = Tsig(sig0,mparams,cparams); 
else
    T_Msig = T_M;
    T_Asig = T_A;
end

if Tini<=T_M
	%No Austenite
	if sig0>sigM
		xplus0=1.0;
	elseif sig0<-sigM
		xminus0=1.0;
	else
		xplus0=0.5;
		xminus0=0.5;
	end
elseif Tini<T_Msig
	%No Austenite
	if sig0>0
		xplus0=1.0;
	elseif sig0<0
		xminus0=1.0;
	else
		xplus0=0.5;
		xminus0=0.5;
	end
elseif Tini<T_Asig
	%Mixed Case
	if Tini==min(Tini,Te) %Heating
		if sig0>0
			xplus0=1.0;
		elseif sig0<0
			xminus=1.0;
		else
			xplus0=0.5;
			xminus0=0.5;
		end
	end
end
%End phase_ini.m