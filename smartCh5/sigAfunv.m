% Jordan E. Massad
% Oct 18, 2004
% Compute Transformation Stresses (Vectorized)
%
%
%SIGAFUNv
%sigA = sigAfunv(T,mparams,cparams)   
%Input 
% T - temperature point/vector
% mparams - material parameters:  
% cparams - computed parameters.
%Output
% sigA - A-M transition stress  
function sigA = sigAfunv(T,mparams,cparams)   

% Initialize Parameters
EM = mparams(2); eT = mparams(3);
dsig = mparams(4);
Er = cparams(1);
du = cparams(2); ds = cparams(3);
Tc = cparams(5);

sigA = T;    %Quick way to make size(sigA)=size(T)
%Don't allow temperatures below second zero (minimum temperature).
T = max(T,cparams(4));

% Nonzero dc Corrections
TR = du/ds;
% dc = 0.00*mparams(6); 
% dc = 0.29*mparams(6);   %Thesis/UCLA
% dc = 0.022*mparams(6);    %UCLA2
% dc = 0.343*mparams(6);  %Ishida
dc = 3e-3*mparams(6); %Ishida NiTiPd
% dc = 0.10*mparams(6);  %Thesis-loops
% dc = 8.0*ds;
Fdc = dc*(T-TR -T.*log(T/TR));

% Reference Chemical Free Energy
dB = du - ds*T + Fdc;
dB = min(dB,EM*eT^2/2);

% Handle T<=Tc (artificially create sigA<0 to calculate sigM=sigA-dsig)
sig1 = (T<=Tc);
y = find(sig1);
sigA(y) = dsig - 2*dB(y)/eT;
% Make sigM reverse direction (decrease) for T<Tc
% a=1.2e-1; sigA(y) = -a*sigA(y);

% Calculate A-M transition stress.
y = find(~(sig1));
A = EM/(1-Er);
sigA(y) = dsig^2/4 + A*(A*eT^2-2*dB(y));
sigA(y) = dsig/2 - A*eT + sqrt(sigA(y)); 

% sigA = dsig^2/4 + EM/(1-Er)*(EM*eT^2/(1-Er)-2*dB);
% sigA = (T>Tc).*(dsig/2 - EM*eT/(1-Er) + sqrt(sigA)); 

%sigA=sigA(:);

%End sigAfun.m