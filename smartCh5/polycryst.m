% Jordan E. Massad
% Oct 18, 2003
% Polycrystal/Homogenization SMA Model
%
%
%POLYCRYST
%[sig,ep,T,xm,xp,tspan] = polycryst(tspan,Tini,mparams,cparams)   
%Input 
% tspan - input time vector
% Tini - initial temperature point
% mparams - material parameters (see matparam.m).
% cparams - computed parameters (see comparam.m).
%Output
% sig - stress input vector
% ep - strain output vector
% T - internal temperature vector
% xm - M- phase fraction vector
% xp - M+ phase fraction vector
% tspan - time vector
function [sig,ep,T,xm,xp,tspan] = polycryst(tspan,Tini,mparams,cparams)

warning off MATLAB:ode15s:IntegrationTolNotMet;
warning off MATLAB:singularMatrix
% Initialize Model Parameters
Teq = cparams(2)/cparams(3); Er_=1-cparams(1); eT=mparams(3); EM=mparams(2);
dsigHi = mparams(13)*(1+1/(Er_/EM/eT*mparams(13)+1));
dsigHi = dsigHi-EM*eT*(mparams(14)-Teq)/Teq/(Er_/EM/eT*mparams(13)+1);
dsigHi = min([EM*eT,dsigHi]);

%Statistical Distribution Parameters
smu = mparams(4);       %Mean Gap Stress [MPa]
svar = mparams(17);     %Gap Variance [MPa^2]
sevar = mparams(18);    %Effective Stress Variance [MPa^2]
%LogNormal Distribution Parameters
lgmu = smu/sqrt(svar/smu^2+1);   %Log Mean
lgvar = log(svar/smu^2+1);       %Log Variance
%Normalization Factors
nfact = 1/sqrt(sevar);           %Effective stress contribution.
nfact = nfact/sqrt(lgvar)/2/pi;  %LogNormal contribution.

%Stdv-based Integration Limits
n = 3;                  %Number of Standard Deviations (n=3 => 99.7% of data)
stdvy = sqrt(sevar);
x0 = smu/exp(sqrt(lgvar))^n; 
x1 = smu*exp(sqrt(lgvar))^n;
x1 = min(dsigHi,x1);    %dsig cannot be too large to guarantee du>0.
disp(['LogNormal integration limits set to (',num2str(x0),',',num2str(x1),').'])
y0=-n*stdvy; y1=n*stdvy;
gw = gweights;
numsub = 5;             %Integration Subintervals - same for both integrations.

%Display Initial Phase Fraction Information
[xplus0, xminus0] = phase_ini(tspan,Tini,mparams,cparams);
disp(['Estimated Initial Phase Fractions:  ','xM+ = ',num2str(xplus0),'  xM- = ',num2str(xminus0)]);

%%Double Integration
tic;
Fy = zeros(4,length(tspan)); FyT=Fy;
Fym = Fy;  Fyp = Fy;
funx = Fy; funy=Fy;    funxT = Fy; funyT=Fy;
funxm = Fy; funym=Fy;  funxp = Fy; funyp=Fy;
ep=zeros(size(tspan)); T=ep; xm=ep; xp=ep;
for jy = 1:numsub
    ysub = gnodes(jy,numsub,y0,y1);
    for inodey = 1:4  %y Subintervals
        sigef = ysub(inodey);
        
        %x Integration
        Fx = repmat(0,size(funx));     FxT = repmat(0,size(funxT));
        Fxm = Fx;     Fxp = Fx;
        for jx = 1:numsub
            xsub = gnodes(jx,numsub,x0,x1);
            for inodex = 1:4 %x Subintervals
                dsig = xsub(inodex);
				[ep,T,xm,xp] = integrand(dsig,sigef,tspan,Tini,mparams,cparams,lgmu,lgvar,sevar);
				funx(inodex,:) = ep; 	funxT(inodex,:) = T;
				funxm(inodex,:) = xm;	funxp(inodex,:) = xp;
			end
            Fx = Fx + funx; %4xlength(tspan) Matrix
            FxT = FxT + funxT; %4xlength(tspan) Matrix
			Fxm = Fxm + funxm;  Fxp = Fxp + funxp; 
        end %End x integration
        
        funy(inodey,:) = gw'*Fx; %1xlength(tspan) Vector
        funyT(inodey,:) = gw'*FxT; %1xlength(tspan) Vector
        funym(inodey,:) = gw'*Fxm;  funyp(inodey,:) = gw'*Fxp; 
    end
    Fy = Fy + funy; %4xlength(tspan) Matrix
    FyT = FyT + funyT; %4xlength(tspan) Matrix
	Fym = Fym + funym;  Fyp = Fyp + funyp; 
end
%Multiply weights and sum:
Fy = (y1-y0)*(x1-x0)/4/numsub^2*gw'*Fy*nfact;
FyT = (y1-y0)*(x1-x0)/4/numsub^2*gw'*FyT*nfact;
Fym = (y1-y0)*(x1-x0)/4/numsub^2*gw'*Fym*nfact;
Fyp = (y1-y0)*(x1-x0)/4/numsub^2*gw'*Fyp*nfact;
disp(['Integration Time:  ',num2str(toc),'s'])
ep = Fy(:);  T = FyT(:);
xm = Fym(:);  xp = Fyp(:);
warning on;

%Re-scale to account for finite interval integration.
refact = Tini/T(1);
T = refact*T; ep = refact*ep;
xm = refact*xm; xp = refact*xp;

%Stress Input
sig = input_sigma(tspan(:));
%External Temperature
Te = input_Te(tspan(:)); 

%Plot Results
plotFun(tspan,sig,Te,ep,T,xm,xp);


%Integrand Function
function [ep,T,xm,xp] = integrand(dsig,sigef,tspan,Tini,mparams,cparams,lgmu,lgvar,sevar) 
mparams(4) = dsig;                   %Update dsig parameter [MPa].
cparams(12) = sigef;                 %Update sigef parameter [MPa].
newcprm = comparam(mparams);         %Update computed parameters.
g = exp((-sigef^2/sevar - log(dsig/lgmu)^2/lgvar)/2)/dsig;
[ep,T,xm,xp] = polyFun(tspan,Tini,mparams,newcprm);
len=length(ep); tlen=length(tspan);
if(len<tlen)
	disp(['ODE15s Vomited after t=',num2str(tspan(len)),'s'])
	disp(['input_sigma=',num2str(input_sigma(tspan(len))),'MPa   ','input_Te=',num2str(input_Te(tspan(len))),'K'])
	disp(['dsig=',num2str(dsig),'MPa   ','sigef=',num2str(sigef),'MPa'])
	disp(['ep=',num2str(ep(end)),'   ','xM=',num2str(xm(end)+xp(end))])
	disp(['T=',num2str(T(len)),'K   ','sigA=',num2str(sigAfun(T(len),mparams,newcprm)),'MPa'])
	disp(['Tmin=',num2str(newcprm(4)),'K   ','T_M=',num2str(newcprm(5)),'K'])
	disp(['Distribution Factor=',num2str(g)])
	disp(' ')
	if(len>0.75*tlen)
		ep(len+1:tlen)=0; T(len+1:tlen)=0; xm(len+1:tlen)=0; xp(len+1:tlen)=0;
	else
		ep=zeros(size(tspan)); T=ep; xm=ep; xp=ep;
	end
end
ep = g*ep;  T = g*T;
xm = g*xm; xp = g*xp;
return;

% End polycryst.m