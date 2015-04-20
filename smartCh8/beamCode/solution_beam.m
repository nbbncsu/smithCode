function bvn = solution_beam(x,N,ell)

h = ell/N;
[N1,N2] = size(x);
[bvn,bvn2] = bevaluate_beam(N,N1,x,h);
