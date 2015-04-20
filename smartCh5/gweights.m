% Jordan E. Massad
% Oct 18, 2004
% 4-pt Gauss-Legendre Quadrature Weights
%
%GWEIGHTS
%gw = gnodes(jsub,numsub,a,b) 
%Input 
% NULL  
%Output
% gw - 4 quadrature weights.  
function gw = gweights

%Standard weights:
w1 = 1/2 + sqrt(30)/36;
w2 = 1/2 - sqrt(30)/36;

gw = [w2;w1;w1;w2];

%  End gweights.m 
