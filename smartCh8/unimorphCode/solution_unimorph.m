function bvn = solution_unimorph(x,N,ell)

h = ell/N;
[N1,N2] = size(x);
[bvn,bvn2] = bevaluate_unimorph(N,N1,x,h);
