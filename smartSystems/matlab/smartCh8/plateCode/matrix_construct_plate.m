function [M,K,C] = matrix_construct_plate(Nx,Nqx,ellx,Ny,Nqy,elly,...
    rho_p,E_p,cD_p,gamma,nu,h,bc);

%  This function constructs the mass, stiffness and damping matrices for the beam.
%  Quadratures are performed using a four point Gaussian rule.

if bc==1
    Nhatnx = Nx+1;
    Nhatny = Ny+1;
    Nhatqx = Nx+1;
    Nhatqy = Ny+1;
else
    Nhatnx = Nx-1;
    Nhatny = Ny-1;
    Nhatqx = Nx-1;
    Nhatqy = Ny-1;
end

hx = ellx/Nx;
hy = elly/Ny; 

hqx = ellx/Nqx;
hqy = elly/Nqy;

gx = gauss_points(Nqx,hqx,0);
gy = gauss_points(Nqy,hqy,0);

% Construct x components

[bvn0,bvn1,bvn2] = bevaluate_plate(Nx,4*Nqx,gx,hx,bc);
[bvq0,bvq1,bvq2] = bevaluate_plate(Nx,4*Nqx,gx,hx,bc);

rpv = gauss_weights(Nqx,hqx);

for n=1:Nhatnx
    bv0n = bvn0(:,n);
    bv1n = bvn1(:,n);
    bv2n = bvn2(:,n);
    for q = 1:Nhatqx
        bv0q = bvq0(:,q);
        bv1q = bvq1(:,q);
        bv2q = bvq2(:,q);
        Mx(q,n) = sum(rpv.*bv0n.*bv0q);
        Kx(q,n) = sum(rpv.*bv2n.*bv2q);
        Wx11(q,n) = sum(rpv.*bv1n.*bv1q);
        Wx02(q,n) = sum(rpv.*bv0n.*bv2q);
        Wx20(q,n) = sum(rpv.*bv2n.*bv0q);
    end
end

% Construct y components

[bvn0,bvn1,bvn2] = bevaluate_plate(Ny,4*Nqy,gy,hy,bc);
[bvq0,bvq1,bvq2] = bevaluate_plate(Ny,4*Nqy,gy,hy,bc);

rpv = gauss_weights(Nqy,hqy);

for n=1:Nhatny
    bv0n = bvn0(:,n);
    bv1n = bvn1(:,n);
    bv2n = bvn2(:,n);
    for q = 1:Nhatqy
        bv0q = bvq0(:,q);
        bv1q = bvq1(:,q);
        bv2q = bvq2(:,q);
        My(q,n) = sum(rpv.*bv0n.*bv0q);
        Ky(q,n) = sum(rpv.*bv2n.*bv2q);
        Wy11(q,n) = sum(rpv.*bv1n.*bv1q);
        Wy02(q,n) = sum(rpv.*bv0n.*bv2q);
        Wy20(q,n) = sum(rpv.*bv2n.*bv0q);
    end
end

D = E_p*(h^3)/(12*(1-nu^2));
cD = cD_p*(h^3)/(12*(1-nu^2));
M = rho_p*h*kron(My,Mx);
K = D*(kron(My,Kx) + kron(Ky,Mx) + 2*kron(Wy11,Wx11) + nu*kron(Wy20,Wx02) ...
       + nu*kron(Wy02,Wx20) - 2*nu*kron(Wy11,Wx11));
C = (cD/D)*K + gamma*kron(My,Mx);
