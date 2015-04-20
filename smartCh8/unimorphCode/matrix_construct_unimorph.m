function [M,K,C] = matrix_construct_unimorph(N,ell,rho,Y,cD,gamma,YI);

%  This function constructs the mass, stiffness and damping matrices for the beam.
%  Quadratures are performed using a four point Gaussian rule.

h = ell/N;
Nhat = N+1;

g = gauss_points(N,h,0);

[bvn,bvn2] = bevaluate_unimorph(N,4*N,g,h);

wt = gauss_weights(N,h);
rpv = wt;

for n=1:Nhat
    bv0n = bvn(:,n);
    bv2n = bvn2(:,n);
    for q = 1:Nhat
        bv0q = bvn(:,q);
        bv2q = bvn2(:,q);
        Spline0(q,n) = sum(rpv.*bv0n.*bv0q);
        Spline2(q,n) = sum(rpv.*bv2n.*bv2q);
    end
end

M = rho*Spline0;
C = gamma*Spline0 + cD*Spline2;
K = YI*Spline2;
