function F = force_plate(Nx,Nqx,ellx,Ny,Nqy,elly,rho_p,E_p,cD_p,gamma,nu,h,bc);

% This function constructs the input vector for the disturbance to the plate.

hx = ellx/Nx;
hy = elly/Ny;

hqx = ellx/Nqx;
hqy = elly/Nqy;

if bc==1
    Nhatnx = Nx+1;
    Nhatny = Ny+1;
else
    Nhatnx = Nx-1;
    Nhatny = Ny-1;
end

gx = gauss_points(Nqx,hqx,0);
gy = gauss_points(Nqy,hqy,0);

% Construct x components

rpv = gauss_weights(Nqx,hqx);

[bvn0,bvn1,bvn2] = bevaluate_plate(Nx,4*Nqx,gx,hx,bc);

for n=1:Nhatnx
    bv0n = bvn0(:,n);
    Fx(n,1) = sum(rpv.*bv0n);
end

% Construct y components

rpv = gauss_weights(Nqy,hqy);

[bvn0,bvn1,bvn2] = bevaluate_plate(Ny,4*Nqy,gy,hy,bc);

for n=1:Nhatny
    bv0n = bvn0(:,n);
    Fy(n,1) = sum(rpv.*bv0n);
end

F = kron(Fy,Fx);
