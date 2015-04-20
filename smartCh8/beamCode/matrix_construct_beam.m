function [M,K,C] = matrix_construct_beam(N,Nq,s,pe1,pe2,ell,...
    rho_b,rho_p,EI_b,EI_p,cD_b,cD_p,gamma_b);

%  This function constructs the mass, stiffness and damping matrices for the beam.
%  Quadratures are performed using a four point Gaussian rule.

h = ell/N;
Nhatn = N+1;
Nhatq = N+1;

for j=1:s+1
    if j==1
        hq = ell/Nq;
        g = gauss_points(Nq,hq,0);
    else
        hq = (pe2-pe1)/Nq;
        g = gauss_points(Nq,hq,pe1);
    end
    [bvn,bvn2] = bevaluate_beam(N,4*Nq,g,h);
    [bvq,bvq2] = bevaluate_beam(N,4*Nq,g,h);

    wt  = gauss_weights(Nq,hq);
    %rpval = ones(size(g));
    %rpv = wt.*rpval;
    rpv = wt;
    
    for n=1:Nhatn
        bv0n = bvn(:,n);
        bv2n = bvn2(:,n);
        for q = 1:Nhatq
            bv0q = bvq(:,q);
            bv2q = bvq2(:,q);
            Spline0(q,n) = sum(rpv.*bv0n.*bv0q);
            Spline2(q,n) = sum(rpv.*bv2n.*bv2q);
        end
    end

    if j==1
        Mb = rho_b*Spline0;
        Kb = EI_b*Spline2;
        Cb = cD_b*Spline2;
        Cg = gamma_b*Spline0;
    else
        Mp = rho_p*Spline0;
        Kp = EI_p*Spline2;
        Cp = cD_p*Spline2;
    end
end

M = Mb + Mp;
K = Kb + Kp;
C = Cb + Cp + Cg;
