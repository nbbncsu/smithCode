function ydot = yprime_beam(t,y,A_mat,F_vec)

force = F_vec*sin(5*2*pi*t);
ydot = A_mat*y + force;
