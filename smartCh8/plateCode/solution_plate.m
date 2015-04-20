function Bvn = solution_plate(gx,Nx,ellx,gy,Ny,elly,bc)

hx = ellx/Nx;
[N1,N2] = size(gx);
[bvnx,bvn1,bvn2] = bevaluate_plate(Nx,N1,gx,hx,bc);

hy = elly/Ny;
[N1,N2] = size(gy);
[bvny,bvn1,bvn2] = bevaluate_plate(Ny,N1,gy,hy,bc);

Bvn = kron(bvny,bvnx);
