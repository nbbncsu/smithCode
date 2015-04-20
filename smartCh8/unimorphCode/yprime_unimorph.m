function ydot = yprime_unimorph(t,y,flag,A_mat,F_vec,B_vec,V)

freq = 1;
Voltage = max(V)*sin(2*freq*pi*t);
switch flag
case ''
    ydot = A_mat*y + Voltage*B_vec;
case 'jacobian'
    ydot = A_mat; 
end




