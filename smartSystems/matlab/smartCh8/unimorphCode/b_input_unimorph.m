function B = b_input_unimorph(N,ell,Kb);

% This function constructs the input from the active component of the unimorph.

h = ell/N;

g = gauss_points(N,h,0);

wt = gauss_weights(N,h);
f1 = wt;

[bvn,bvn2] = bevaluate_unimorph(N,4*N,g,h);

for n=1:N+1
    bv2n = bvn2(:,n);
    B(n,1) = Kb*sum(f1.*bv2n);
end
