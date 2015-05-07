function F = force_unimorph(N,ell,fval);

% This function constructs the input vector for the disturbance to the
% beam.  The force is currently coded for a point force at the end of
% the beam to be used for computing blocked forces.

h = ell/N;

g = gauss_points(N,h,0);

wt = gauss_weights(N,h);
f1 = wt;

[bvn,bvn2] = bevaluate_unimorph(N,4*N,g,h);

for n=1:N+1
    bv0n = bvn(:,n);
    %F(n,1) = sum(f1.*bv0n);
end

[bvq,bvq2] = bevaluate_unimorph(N,1,ell,h);
F = fval*bvq';
