% Jordan E. Massad
% Oct 18, 2004
% Composite 4-pt Gauss-Legendre Quadrature Nodes
% Equilength Subintervals on [a,b]
%
%GNODES
%xj = gnodes(jsub,numsub,a,b) 
%Input 
% jsub - index of jth subinterval.
% numsub - number of subintervals.
% a - start of main interval  
% b - end of main interval  
%Output
% xj - 4 quadrature nodes in subinterval jsub  
function xj = gnodes(jsub,numsub,a,b)

if jsub<=0
    error('Nonpositive subinterval index.')
end
if jsub>numsub
    error('Subinterval index exceeds number of subintervals.')
end
%Standard nodes on [-1,1]:
x1 = sqrt(15-2*sqrt(30))/sqrt(35);
x2 = sqrt(15+2*sqrt(30))/sqrt(35);
x = [-x2;-x1;x1;x2];

%Nodes mapped to subinterval jsub:
xj = a + (b-a)/2/numsub.*(2*jsub-1 + x);

%  End gnodes.m 
