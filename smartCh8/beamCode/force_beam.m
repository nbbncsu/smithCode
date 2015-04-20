function F = force_beam(N,Nq,ell);

% This function constructs the input vector for the disturbance to the
% beam.  It is currently coded for a uniform (in space) force such as
% might be generated by a pressure field.

h = ell/N;
hq = ell/Nq;

g = gauss_points(Nq,hq,0);
f1 = gauss_weights(Nq,hq);

[bvn,bvn2] = bevaluate_beam(N,4*Nq,g,h);

for n=1:N+1
    bv0n = bvn(:,n);
    F(n,1) = sum(f1.*bv0n);
end
