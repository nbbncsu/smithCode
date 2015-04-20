% Jordan E. Massad
% Oct 18, 2004
% Plotting Function 
%
%
%PLOTFUN
%[] = plotFun(t,ep,T,xm,xp)   
%Input 
% t - input time vector
% sig - input stress vector
% Te - input external temperature vector
% ep - strain output vector
% T - internal temperature vector
% xm - M- phase fraction vector
% xp - M+ phase fraction vector
function [] = plotFun(t,sig,Te,ep,T,xm,xp)

figure(1);
if (norm(diff(sig))>0 & norm(diff(Te))>0)
    subplot(2,1,1), plot(100*ep,sig,'.-')
    ylabel('Stress (MPa)')
	xlabel('Strain (%)')
    subplot(2,1,2), plot(Te,100*ep,'.-')
	ylabel('Strain (%)')
	xlabel('External Temperature (K)')
elseif norm(diff(sig))>0
    plot(100*ep,sig,'.-')
    ylabel('Stress (MPa)')
	xlabel('Strain (%)')
elseif norm(diff(Te))>0
    plot(Te,100*ep,'.-')
	ylabel('Strain (%)')
	xlabel('External Temperature (K)')
else
	plot(T,100*ep,'.-')
	ylabel('Strain (%)')
	xlabel('Internal Temperature (K)')
end
title('Hysteresis')

figure(2);
subplot(2,1,1)
%plot(1e-6*sig,epA,'r',1e-6*sig,epMp,'b',1e-6*sig,epMm,'g')
plot(t,xm,'b:',t,xp,'b',t,1-xm-xp,'r')
axis([0 t(end) -0.05 1.05]);
xlabel('Time (s)')
ylabel('Phase Fraction')
subplot(2,1,2)
plot(t,T)
axis tight;
xlabel('Time (s)')
ylabel('Temperature (K)')
hold off;

% End plotFun.m