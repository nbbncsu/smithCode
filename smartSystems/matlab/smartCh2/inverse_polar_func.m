function [E, P_out, Ec_int, Ei_int, Ec_gpts, Ei_gpts, nu_one, nu_two] = ...
    inverse_polar_func(params, P, Delta_t, freq_dep_yes, density_type)

%%% Set the number of polarization points and the number of quadrature intervals.

Nk = length(P);
Np = 20;
Nq = 20;
Ni = 4*Np;
Nj = 4*Nq;
   
%%% Extract the parameters.
%%% In the case of general densities, the integration limits are fixed
%%% parameters. Otherwise they are determined dynamically from the
%%% lognormal/normal or normal/normal parameters.

if density_type == 3
    Ec_int = params(1);
    Ei_int = params(2);
    params = params(3:end);
end;

eta = params(1);
P_R = 1;
if density_type == 1
    C = params(2);
    b_one = params(3);
    c_one = sqrt(C);
    b_two = params(4);
    c_two = sqrt(C);
    Ec_int = sqrt(b_one*7);
    Ei_int = sqrt(b_two*7);
elseif density_type == 2
    C = params(2);
    Ec_bar = params(3);
    c = params(4);
    c_one = sqrt(C);
    b_two = params(5);
    c_two = sqrt(C);
    Ec_int = Ec_bar*200^c;
    Ei_int = sqrt(b_two*7);
end
           
%%%  Setup the quadrature points Ec_gpts and Ee_gpts along
%%%  with the densities nu_one and nu_two.

hp = Ec_int/Np;
v = gauss_weights(Np,hp);
Ec_gpts = gauss_points(Np,hp,0);
if density_type == 1
    nu_one = c_one*exp(-((Ec_gpts).^2)/b_one);    
elseif density_type == 2
    nu_one = c_one*exp(-(log(Ec_gpts/Ec_bar)/(2*c)).^2);
elseif density_type == 3
    nu_one = params(2:2+Ni-1);
end

hq = 2*Ei_int/Nq; % factor of two b/c this is neg-inf to pos-inf integrand
w = gauss_weights(Nq,hq);
Ei_gpts = gauss_points(Nq,hq,0-Ei_int);
if density_type == 1 | density_type == 2 | density_type == 4
    nu_two_pos = c_two*exp(-((Ei_gpts(Nj/2+1:Nj)).^2)/b_two);
elseif density_type == 3
    nu_two_pos = params(2+Ni:2+Ni+Nj/2-1);
end
nu_two_neg = rot90(rot90(nu_two_pos));
nu_two = [nu_two_neg; nu_two_pos];

%%% Setup the initial vectors V and W and the matrices epsilon_c and epsilon_k.

V = v.*nu_one;
W = w.*nu_two;
eps_c = Ec_gpts*ones(1,Nj);
eps_k = ones(Ni,1)*Ei_gpts';

%%% Determine the electric field and the output polarization.

E = zeros(Nk,1);
P_out = zeros(Nk,1);

if freq_dep_yes == 0
    Delta = ones(Ni,Nj);
    Delta(:,1:Nj/2) = -1;

    E(1) = 0;

    step_init = 1.6e6; % intial step size
    step(1) = step_init;
    sct_prev = 1; % sct and sct_prev count the number of steps on a given loop

    for k=2:Nk

        %%% Set the step size adaptively based on the number
        %%% of steps taken on the previous loop.

        if sct_prev == 1
            step(k) = step_init;
        elseif sct_prev >= 10 
            step(k) = 2*step(k-1);
        else 
            step(k) = step(k-1);
        end;

        dP = P(k) - P(k-1);
        dE = step(k)*dP;

        E_temp = E(k-1);
        P_temp(1) = P(k-1);
        Delta_temp = Delta;
        eps_k_temp = eps_k;

        sct = 1;

        %%% If there is a change in P, step through changes in E
        %%% until the desired P value is passed. Then interpolate
        %%% to get an estimate of the correct E value.
        
        if dP == 0;
            E(k) = E(k-1);
        else
            while sign(dP)*(P(k) - P_temp(sct)) >= 0;
                E_temp = E_temp + dE;
                eps_k_temp = eps_k_temp + dE;
                Delta_temp = sign(eps_k_temp + eps_c.*Delta_temp);
                P_bar = (1/eta)*eps_k_temp + P_R*Delta_temp;
                P_temp(sct+1) = V'*P_bar*W;
                sct = sct+1;
            end;
            slope = (P_temp(sct) - P_temp(sct-1))/dE;
            E(k) = E_temp - (P_temp(sct) - P(k))/slope;
            sct_prev = sct;
        end;

        %%% Determine the output polarization value.
        
        dE = E(k) - E(k-1);
        eps_k = eps_k + dE;
        Delta = sign(eps_k + eps_c.*Delta);
        P_bar = (1/eta)*eps_k + P_R*Delta;
        P_out(k) = V'*P_bar*W;
    end;
end;
