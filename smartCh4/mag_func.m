function [P, Ec_int, Ei_int, Ec_gpts, Ei_gpts, nu_one, nu_two] = ...
    polar_func(params, E, Delta_t, freq_dep_yes, density_type)

%%% Set the number of field points and the number of quadrature intervals.

Nk = length(E);
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

if freq_dep_yes == 1
    jump_freq = params(1); % jump_freq = 1/Tau
    beta_param = params(2); % beta_param = sqrt(eta*V/k*T)
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
if density_type == 1 | density_type == 2
    nu_two_pos = c_two*exp(-((Ei_gpts(Nj/2+1:Nj)).^2)/b_two);
elseif density_type == 3
    nu_two_pos = params(2+Ni:2+Ni+Nj/2-1);
end
nu_two_neg = rot90(rot90(nu_two_pos));
nu_two = [nu_two_neg; nu_two_pos];

%%% Setup the initial vectors V and W and the matrices epsilon_c and epsilon_k

V = v.*nu_one;
W = w.*nu_two;
eps_c = Ec_gpts*ones(1,Nj);
eps_k = ones(Ni,1)*Ei_gpts';

%%% Determine the polarization.

P = zeros(Nk,1);

if freq_dep_yes == 0
    Delta = ones(Ni,Nj);
    Delta(:,1:Nj/2) = -1;

    P_bar = (1/eta)*eps_k + P_R*Delta;
    P(1) = V'*P_bar*W;

    %%%  Start iteration

    for k=2:Nk
        dE = E(k) - E(k-1);
        eps_k = eps_k + dE;
     
        Delta = sign(eps_k + eps_c.*Delta);

        P_bar = (1/eta)*eps_k + P_R*Delta;
        P(k) = V'*P_bar*W; 
    end;
elseif freq_dep_yes == 1
    P_I = P_R - eps_c/eta;
    x_pos_prev = ones(Ni,Nj);
    x_pos_prev(:,1:Nj/2) = 0;
    epsilon = 1e-3;
    fudge = 1e-50;

    %  Start iteration

    for k = 2:Nk
        E_cur = E(k) + eps_k;
        P_cur_min_L = E_cur/eta - P_R;
        P_cur_min_R = E_cur/eta + P_R;
        r_neg_pos = (sign(eps_c - E_cur) + 1)/2.*...
            (erfc((P_cur_min_L - (-P_I + epsilon))/(beta_param*sqrt(2))) -...
            erfc((P_cur_min_L - (-P_I - epsilon))/(beta_param*sqrt(2))))./...
            (erfc((P_cur_min_L - (-P_I + epsilon))/(beta_param*sqrt(2))) + fudge) + (1 - sign(eps_c - E_cur))/2;
        r_pos_neg = (sign(eps_c + E_cur) + 1)/2.*...
            (erfc((P_I - epsilon - P_cur_min_R)/(beta_param*sqrt(2))) -...
            erfc((P_I + epsilon - P_cur_min_R)/(beta_param*sqrt(2))))./...
            (erfc((P_I - epsilon - P_cur_min_R)/(beta_param*sqrt(2))) + fudge) + (1 - sign(eps_c + E_cur))/2;
        p_neg_pos = jump_freq*r_neg_pos;
        p_pos_neg = jump_freq*r_pos_neg;
        x_pos_cur = (x_pos_prev + Delta_t*p_neg_pos)./(1 + Delta_t*(p_pos_neg + p_neg_pos));
        x_pos_prev = x_pos_cur;
        P_neg_avg = -beta_param/sqrt(2*pi)*exp(-(-P_I - epsilon - P_cur_min_L).^2/(2*beta_param^2))./...
            (0.5*erfc(-(-P_I - epsilon - P_cur_min_L)/(beta_param*sqrt(2))) + fudge) + P_cur_min_L;
        P_pos_avg = beta_param/sqrt(2*pi)*exp(-(P_I + epsilon -P_cur_min_R).^2/(2*beta_param^2))./...
            (0.5*erfc((P_I + epsilon - P_cur_min_R)/(beta_param*sqrt(2))) + fudge) + P_cur_min_R;
        P_bar = x_pos_cur.*P_pos_avg + (1 - x_pos_cur).*P_neg_avg;
        P(k) = V'*P_bar*W;
    end;
end;
