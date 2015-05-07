% Jordan E. Massad
% Oct 18, 2004
% Macroscopic Strain Averages
% Vectorized Routine (It's Fast!!)
%   2003 Changes 1)T<Tc behavior 2)No sig>sigM|sigA. 3)Use erfcx. 4)Fix epA. 
%
%MEPSv
%[epA,epMplus,epMminus] = mepsv(sig,T,mparams,cparams) 
%Input 
% sig - stress vector
% T - temperature vector
% mparams - material parameters:  
% cparams - computed parameters.
%Output
% epA,epMplus,epMminus - average strain vectors
function [epA,epMplus,epMminus] = mepsv(sig,T,mparams,cparams)    

% Initialize Parameters
EA = mparams(1);
EM = mparams(2);
eT = mparams(3);
Tc = cparams(5);
epA=zeros(size(sig)); epMminus=epA; epMplus=epA; 
expo=epA; pAm=epA; pAp=epA;

%Don't allow temperatures below Helmholtz minimum temperature.
T = max(T,cparams(4));

qAT = cparams(8)./sqrt(T);
qMT = cparams(9)./sqrt(T);

%Transition stress.
sigA = sigAfunv(T,mparams,cparams);
sigM = sigA-mparams(4);

%Standard Deviations of Numerators
wA=1/sqrt(2)./qAT;  wM=1/sqrt(2)./qMT; 
z = 10;   %Number of standard deviations for pwise approx.

%<eA>
%if (abs(sig)>=sigA+z/2*wA)
sig1 = (abs(sig)>=sigA).*(T>Tc);
y = find(sig1);
expo(y) = 2*qAT(y).*sigA(y); 
epA(y) = sign(sig(y)).*(exp(-expo(y).^2)-1)./erf(expo(y))./qAT(y)/EA/sqrt(pi);
epA(y) = epA(y) + sign(sig(y)).*sigA(y)/EA;
% epA(y) = sign(sig(y)).*sigA(y)/EA;
%elseif (abs(sig)<=sigA-z*wA)
sig2 = ~sig1.*(abs(sig)<=sigA-z*wA).*(T>Tc);
y = find(sig2);
epA(y) = sig(y)/EA;
%else
y = find(~(sig1+sig2).*(T>Tc));
expo(y) = qAT(y).*(sigA(y)+sig(y)); 
% pAm(y) = 1./(erfcx(qAT(y).*(sig(y)-sigA(y))).*exp(4*qAT(y).^2.*sigA(y).*sig(y))-erfcx(expo(y)));
% expo(y) = qAT(y).*(sigA(y)-sig(y)); 
% pAp(y) = 1./(erfcx(-expo(y))-erfcx(qAT(y).*(sigA(y)+sig(y))).*exp(-4*qAT(y).^2.*sigA(y).*sig(y)));
% epA(y) = wA(y)/EA*sqrt(2/pi).*(pAm(y) - pAp(y));
depA = (erf(qAT(y).*(sigA(y)-sig(y))) + erf(expo(y))).*(EA*qAT(y)*sqrt(pi));
epA(y) = (exp(-expo(y).^2) - exp(-(qAT(y).*(sig(y)-sigA(y))).^2));
epA(y) = epA(y)./depA;
epA(y) = epA(y) + sig(y)/EA;
%end

%<eMplus>
% if (sig<=sigM-z*wM)
sig1 = (sig<=sigM-z*wM);
y = find(sig1);
epMplus(y) = eT+sigM(y)/EM;% +sqrt(2/pi)/EM*wM(y);
% elseif (sig<sigM+z/2*wM)
sig2 = ~sig1.*(sig<sigM);
y = find(sig2);
expo(y) = qMT(y).*(sigM(y)-sig(y));
epMplus(y) =  wM(y)/EM*sqrt(2/pi)./erfcx(expo(y));
epMplus(y) = epMplus(y) + sig(y)/EM + eT;
% depM = erfc(expo(y)).*(EM*qMT(y)*sqrt(pi));
% epMplus(y) = (exp(expo(y).^2).*depM); 
% epMplus(y) = 1./epMplus(y) + sig(y)/EM + eT;
% else
y = find(~(sig1+sig2));
epMplus(y) = eT+sig(y)/EM; 
% end

%<eMminus>
% if (sig<=-sigM-z/2*wM)
sig1 = (sig<=-sigM);
y = find(sig1);
epMminus(y) = -eT+sig(y)/EM;
% elseif (sig<-sigM+z*wM)
sig2 = ~sig1.*(sig<-sigM+z*wM);
y = find(sig2);
expo(y) = qMT(y).*(sig(y)+sigM(y));
epMminus(y) =  -wM(y)/EM*sqrt(2/pi)./erfcx(expo(y));
epMminus(y) = epMminus(y) + sig(y)/EM - eT;
% depM = erfc(qMT(y).*(sigM(y)+sig(y))).*(EM*qMT(y)*sqrt(pi));
% epMminus(y) = (-exp((qMT(y).*(sig(y)+sigM(y))).^2).*depM); 
% epMminus(y) = 1./epMminus(y) + sig(y)/EM - eT;   
% else
y = find(~(sig1+sig2));
epMminus(y) = -eT-sigM(y)/EM;% -sqrt(2/pi)/EM*wM(y);
% end

%End mepsv.m