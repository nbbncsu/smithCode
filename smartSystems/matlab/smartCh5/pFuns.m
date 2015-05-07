% Jordan E. Massad
% Oct 18, 2004
% Compute Transformation Likelihoods (all temperatures)
%
%
%PFUNS
%[Pm,Pp,pAm,pAp] = pFuns(sig,T,mparams,cparams)   
%Input 
% sig - input stress point.
% T - temperature point.
% mparams - material parameters.  
% cparams - computed parameters.  
%Output
% Pm,Pp,pAm,pAp - likelihoods
function [Pm,Pp,pAm,pAp] = pFuns(sig,T,mparams,cparams)   

% Initialize Parameters
T_A = cparams(5);
sqEr = sqrt(cparams(1));
z = 5;   %Number of standard deviations for pwise approx.
pAm = 1; pAp = 1;

%Don't allow temperatures below Helmholtz minimum temperature.
%%Bad ODE iterate gets unrealistic, low temperatures.
T = max(T,cparams(4));

%Exponent Factors
qAT = cparams(8)/sqrt(T);
qMT = cparams(9)/sqrt(T);

%Standard Deviations
wA=1/sqrt(2)/qAT;  wM=1/sqrt(2)/qMT;

%Transition stress.
sigA = sigAfun(T,mparams,cparams);
sigM = sigA-mparams(4);

if ((T>T_A)&(sigA>1e-2))
	%pA-
	if sig<-sigA
		pAm = min(1e4,1/erf(2*qAT*sigA));
	elseif sig<=(z*wA-sigA)
		expo = qAT*(sigA+sig);
		pAm = 1/(erfcx(qAT*(sig-sigA))*exp(4*qAT^2*sigA*sig)-erfcx(expo));
	else
		pAm = 0;
	end
	
	%pA+
	if sig<(sigA-z*wA)
		pAp = 0;
	elseif sig<=sigA
		expo = qAT*(sigA-sig);
		pAp = 1/(erfcx(-expo)-erfcx(qAT*(sigA+sig))*exp(-4*qAT^2*sigA*sig));
	else
		pAp = min(1e0,1/erf(2*qAT*sigA));
	end
end

%P-
if sig<(-z*wM-sigM)
	Pm = 0;
elseif sig<=-sigM
	expo = qMT*(sigM+sig);
	Pm =  sqEr/erfcx(expo);
else
    Pm = sqEr;
end

%P+
if sig<sigM
    Pp = sqEr;
elseif sig<=(z*wM+sigM)
	expo = qMT*(sigM-sig);
	Pp =  sqEr/erfcx(expo);
else
    Pp = 0;
end

% End pFuns.m