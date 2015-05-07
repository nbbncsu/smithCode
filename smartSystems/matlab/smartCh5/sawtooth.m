% Jordan E. Massad
% Oct 18, 2004
% Sawtooth Stress/Temperature Input
%
%
%SAWTOOTH
%sigma = sawtooth(peak,cycper,t)   
%Input 
% peak - cycle amplitude.
% cycper - cycle period.
% t - time value/vector.
%Output
% sigma - stress value/vector.
function sigma = sawtooth(peak,cycper,t)

rate = peak/cycper*2;  %[MPa/s] or [K/s]

%%%Linear interpolation implementation
si = [0 1 0]'*peak;
ti = [0 cumsum(abs(diff(si')))/rate]';
t = mod(t,ti(end));
sigma = interp1(ti,si,t);

%End sawtooth.m