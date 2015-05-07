function ydot = yprime_plate(t,y,flag,A_mat,F_vec,freq_set)

switch flag
case ''
    force = F_vec*sum(sin(freq_set*2*pi*t));
    ydot = A_mat*y + force;

case 'jacobian'
    ydot = A_mat;
end

