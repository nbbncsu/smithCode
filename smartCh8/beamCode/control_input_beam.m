function B = control_input_beam(N,Nq,ell,pe1,pe2,Kb);

h = ell/N;
hq = (pe2-pe1)/Nq;

g = gauss_points(Nq,hq,pe1);
f1 = gauss_weights(Nq,hq)*Kb;

[bvn,bvn2] = bevaluate_beam(N,4*Nq,g,h);

for n=1:N+1
    bv2n = bvn2(:,n);
    B(n,1) = sum(f1.*bv2n);
end
