function [ellx,elly,actuator,act_sense,aval,kval,Ms,lams,alpha,cval,rho_p,...
    E_p,cD_p,gamma,nu,h,Kb,Nx,Ny,Nqx,Nqy,N_act,N_newt,N,bc] = parameters_plate;

Nx = 4;
Ny = Nx;
Nqx = Nx;
Nqy = Nx;
N_act = 4;
N_newt = 16;
bc = 1;
if bc==1
  N = (Nx+1)*(Ny+1);
else
  N = (Nx-1)*(Ny-1);
end

h = .0016;              % (m)
ellx = 0.4;
elly = 0.6;

rho_p = 2700;                 
E_p = 4.1e+10;                  
cD_p = 2.5e+5;              
gamma = .18;            % (s N/m^2)
nu = .33;
ell_r = .0254;          % (m)     
          
if bc==1                                            %  ------------------
    p1 = [0.10,  0.11,  0.395, 0.405];              %  |                | 
    p2 = [0.10,  0.11,  0.195, 0.205];              %  | p1             |
    p3 = [0.145, 0.155, 0.10,  0.11 ];              %  |                |
    p4 = [0.245, 0.255, 0.10,  0.11 ];              %  | p2             |
    actuator = [p1; p2; p3; p4];                    %  |     p3   p4    |
    act_sense = [1 1 0 0                            %  ------------------
                 0 0 1 1];                              
else
    p1 = [0.10,  0.11,  0.295, 0.305];              %  ------------------  
    p2 = [0.195, 0.205, 0.10,  0.11];               %  |       p4       |
    p3 = [0.29,  0.30,  0.295, 0.305];              %  |                |  
    p4 = [0.195, 0.205, 0.49,  0.50];               %  | p1          p3 |
    actuator = [p1; p2; p3; p4];                    %  |                |
    act_sense = [1 0 -1 0                           %  |       p2       |
                 0 1 0 -1];                         %  ------------------
end

E_H = 7.0e+10;                                      
cD_H = 8.0e+5;
A_mag = .0064;

rho_r = .433;
EI_r = 0.0;
cD_r = 0.0;

Ms = 1.3236e+5;
aval = 7105;
kval = 7002;
alpha = .007781;
cval = .3931;
lams = (2/5)*90e-6 + (3/5)*1600e-6;

Kb = 3*(lams/Ms^2)*A_mag*E_H*(h/2 + ell_r)^2;
