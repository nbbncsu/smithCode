function [bv,bv2] = bevaluate_unimorph(N,Np,X,h)

h3 = 1/(h^3);
bval = zeros(Np,N+3);
bval2 = zeros(Np,N+3);
one = 1;
for ell=1:Np
    x = X(ell,1);
    for kk = 1:N+3
        k = kk-2;
        ind1 = k-2;
        ind2 = k-1;
        ind3 = k;
        ind4 = k+1;
        ind5 = k+2;
        if (ind1*h<=x)&(x<ind2*h)
            b1 = h3*((x-ind1*h*one).^3);
            bval(ell,kk) = b1;
            b1 = h3*(6*(x-ind1*h*one));
            bval2(ell,kk) = b1;
        elseif (ind2*h<=x)&(x<ind3*h)
            b2 = h3*(h^3*one + 3*h^2*(x-ind2*h*one) + 3*h*(x-ind2*h*one).^2);
            b2 = b2 - h3*3*((x-ind2*h*one).^3);
            bval(ell,kk) = b2;
            b2 = h3*(6*h*one - 18*(x-ind2*h*one));
            bval2(ell,kk) = b2;
        elseif (ind3*h<=x)&(x<ind4*h) 
            b3 = h3*(h^3*one + 3*h^2*(ind4*h*one-x) + 3*h*(ind4*h*one-x).^2);
            b3 = b3 - h3*3*((ind4*h*one-x).^3);
            bval(ell,kk) = b3;
            b3 = h3*(6*h*one - 18*(ind4*h*one-x));
            bval2(ell,kk) = b3;
        elseif (ind4*h<=x)&(x<ind5*h)
            b4 = h3*((ind5*h*one-x).^3);
            bval(ell,kk) = b4;
            b4 = h3*(6*(ind5*h*one-x));
            bval2(ell,kk) = b4;
        else
            bval(ell,kk) = 0;
            bval2(ell,kk) = 0;
        end
    end
end

bv = zeros(Np,N+1);
bv2 = zeros(Np,N+1);
bv(:,1) = bval(:,2) - 2*bval(:,1) - 2*bval(:,3);
bv(:,2:N+1) = bval(:,4:N+3);
bv2(:,1) = bval2(:,2) - 2*bval2(:,1) - 2*bval2(:,3);
bv2(:,2:N+1) = bval2(:,4:N+3);
