% Jordan E. Massad
% Oct 18, 2004
% Stress Input
%
%
%INPUT_SIGMA
%sigma = input_sigma(t)   
%Input 
% t - time value/vector.
%Output
% sigma - stress value/vector.
function sigma = input_sigma(t)   

rate = 4/3;             %Stress Rate [MPa/s]
sigmax = 200;           %Stress Amplitide [MPa]
cycper = 2*sigmax/rate; %Cycle Period [s]

%Sinusoidal
sigma = sigmax*sin(2*pi*t/cycper/2);

%SawTooth (piecewise linear)
% sigma = sawtooth(sigmax,cycper,t);

%Constant
% sigma = repmat(200e0,size(t));


%  End input_sigma.m