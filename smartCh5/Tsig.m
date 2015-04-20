% Jordan E. Massad
% Oct 18, 2004
% Compute Stress-dependent Transformation Temperatures
%
%
%TSIG
%[Tminsig, T_Msig, T_Asig, ] = Tsig(sig,mparams,cparams)   
%Input 
% mparams - material parameters
% cparams - computed parameters
%Output
% T_Msig - A-M transformation temperature  
% T_Asig - M-A transformation temperature
function [T_Msig, T_Asig] = Tsig(sig,mparams,cparams)    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eT=mparams(3); dsig=mparams(4); 
T_M=mparams(11); T_A=mparams(12);
s_lo=mparams(15); s_hi=mparams(13);
T_lo=mparams(16); T_hi=mparams(14);
Er=cparams(1); EM=mparams(2);
Teq = mparams(12);
dc = mparams(19); 

a = dc/cparams(3);   %Nondimensional Chemical Free Energy Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Stress-dependent Temperatures (nonzero dc)
opts=optimset('fzero'); opts=optimset(opts,'disp','off');

%%T_M [K] (sigA=0 2-well to 3-well temperature)
dBT_M = dsig*eT/2 - F(sig,eT,Er,EM,dsig);
T_Msig = Teq*fzero(@dBndim,1.0,opts,cparams(2),a,dBT_M);

%%T_A [K] (wM=eT superelastic boundary)
dBT_A = dsig*eT/2 - F(sig+dsig,eT,Er,EM,dsig);
T_Asig = Teq*fzero(@dBndim,1.0,opts,cparams(2),a,dBT_A);

%%T_Asup [K] (wA=eT high temp beyond superelasticity)
% dB3=eT/2/cparams(1)*(dsig - ((mparams(1)+mparams(2))*eT));
% T_Asup = fzero(@dBfun,Teq,opts,cparams(2),cparams(3),dc,Teq,dB3);
% Tmin = Teq*fzero(@dBndim,[eps exp(-1/a)],opts,cparams(2),a,0)

% End Tsig()

function Fsig = F(sig,eT,Er,EM,dsig)
Fsig = sig*(eT + (1-Er)/2*(sig-dsig)/EM);

%%Chemical Free Energy Difference
function dB = dBfun(T,du,dh,dc,TR,E)
dB = du -T*dh +dc*(T-TR-T.*log(T/TR)) - E;

%%Chemical Free Energy Difference (nondimensional)
function db = dBndim(theta,du,a,E)
db = 1 - theta + a*(theta-1 - theta.*log(theta)) - E/du;

%  End of Tsig.m