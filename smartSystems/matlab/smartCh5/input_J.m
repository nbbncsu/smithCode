% Jordan E. Massad
% Oct 18, 2004
% External Joule Heating Input
%
%
%INPUT_J
%QJ = input_J(t)   
%Input 
% t - time value/vector.
%Output
% QJ - Joule heating density rate [W/m^3].
function Qj = input_J(t)   

% I = 28e-1;  %Constant Amperage [A]
% Xa = 8e-6*0.1778e-2; %Cross-sectional area [m^2]
% R = (I/Xa)^2;   %Joule heating factor.

V = 5.0;        %Constant Voltage [V] 
L = 1.496e-2;   %Conductive Path Length [m]
R = (V/L)^2;    %Joule Heating Factor

% twait - secs before applying Joule source
% ton - turn on source time. 
% toff - shut off source time. 
% Repeat periodically until tfinal.

% ton=0.021; toff=0.031; tfinal=0.11; %2V inner loops
% ton=0.020; toff=0.023; tfinal=0.023; %5V fast
% t = mod(t,toff).*(t<tfinal);  
% Qj = repmat(R,size(t)).*(t>ton);

%%No Joule Input
Qj = 0;

%End input_J.m