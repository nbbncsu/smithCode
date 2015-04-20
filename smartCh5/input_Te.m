% Jordan E. Massad
% Oct 18, 2004
% External Temperature Input
%
%
%INPUT_TE
%Te = input_Te(t)   
%Input 
% t - time value/vector.
% Tini - initial environment temperature [K].
%Output
% Te - external temperature value/vector [K].
function Te = input_Te(t)   

Tini = 403;         %Initial Temperature [K] 
Tend = 203;         %Final Temperature [K] 
Tmax = Tend-Tini;   %Heating/cooling range [K].
hper = Tmax*12/8;   %Heating Period 10/60 [s]

% Te = repmat(Tini,size(t));
% Te = Tini + 10*(t<20);
% Te = 100*(sin(2*pi/90*t)).*(t<23) + Tini;
% r = 0.17;  %Heating/Coling Rate [K/s]
% tcycle = abs(Tend-Tini)/r;
% t = mod(t,2*tcycle);
% tcool = find(t<=tcycle);
% theat = find(t>tcycle);
% Te(tcool) = Tini - r*t(tcool);
% Te(theat) = 2*Tend-Tini + r*t(theat);

%SawTooth (piecewise linear)
Te = sawtooth(Tmax,hper,t) + Tini;

%End input_Te.m