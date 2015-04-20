%
%         afm_model_demo.m
%

clear all
close all
format short e

%%% Set options.

freq_dep_yes = 0; % must be 0 or 1, 1 for thermal activation
density_type = 1; % 1 for normal/normal, 3 for general densities

%%% Load the model parameters.

if freq_dep_yes == 0
    if density_type == 1
        load params_afm_nor_demo.dat
        params = params_afm_nor_demo;
    elseif density_type == 3
        load params_afm_gen_demo.dat
        params = params_afm_gen_demo;
    end;
elseif freq_dep_yes == 1
    if density_type == 1
        load params_afm_nor_freq_demo.dat
        params = params_afm_nor_freq_demo;
    elseif density_type == 3
        load params_afm_gen_freq_demo.dat
        params = params_afm_gen_freq_demo;
    end;
end;

%%% Set the value of Delta_t (which determines the thermal activation level).

Delta_t = 1e-3;
freq = 1/(140*Delta_t);

%%% Setup the input field.

field_1a = linspace(0,7000,71)';
field_1b = linspace(7000,0,71)';
field_1 = [field_1a(1:end-1); field_1b];
field_2a = linspace(500,6500,61)';
field_2b = linspace(6500,500,61)';
field_2 = [field_2a(1:end-1); field_2b];
field_3a = linspace(1000,6000,51)';
field_3b = linspace(6000,1000,51)';
field_3 = [field_3a(1:end-1); field_3b];
field_4a = linspace(1500,5500,41)';
field_4b = linspace(5500,1500,41)';
field_4 = [field_4a(1:end-1); field_4b];
field_5a = linspace(2000,5000,31)';
field_5b = linspace(5000,2000,31)';
field_5 = [field_5a(1:end-1); field_5b];
field_6a = linspace(2500,4500,21)';
field_6b = linspace(4500,2500,21)';
field_6 = [field_6a(1:end-1); field_6b];
field_7a = linspace(3000,4000,11)';
field_7b = linspace(4000,3000,11)';
field_7 = [field_7a(1:end-1); field_7b];

bridge_1 = linspace(0,500,6)';
bridge_2 = linspace(500,1000,6)';
bridge_3 = linspace(1000,1500,6)';
bridge_4 = linspace(1500,2000,6)';
bridge_5 = linspace(2000,2500,6)';
bridge_6 = linspace(2500,3000,6)';

field = [field_1(1:end-1); field_1(1:end-1); bridge_1(1:end-1);...
         field_2(1:end-1); field_2(1:end-1); bridge_2(1:end-1);...
         field_3(1:end-1); field_3(1:end-1); bridge_3(1:end-1);...
         field_4(1:end-1); field_4(1:end-1); bridge_4(1:end-1);...
         field_5(1:end-1); field_5(1:end-1); bridge_5(1:end-1);...
         field_6(1:end-1); field_6(1:end-1); bridge_6(1:end-1);...
         field_7(1:end-1); field_7];

clear field_1a field_1b field_2a field_2b field_3a field_3b field_4a field_4b
clear field_5a field_5b field_6a field_6b field_7a field_7b
clear bridge_1 bridge_2 bridge_3 bridge_4 bridge_5 bridge_6 bridge_7
 
%%% Run the polarization model. The two values prepended to the field saturate
%%% the model. The resulting values of P are then removed.

[P, Ec_int, Ei_int, Ec_gpts, Ei_gpts, nu_one, nu_two] = ...
    polar_func(params(4:end), [0; 1e7; field], Delta_t,...
        freq_dep_yes, density_type);
P = P(3:end);

%%% Run the spring model. The disp_shift parameter is needed because
%%% the output measurement device is not calibrated to have a 0 reading
%%% correspond to 0 displacement.

disp_shift = params(1);
k = params(2);
c = params(3);
w = k;
A=[0 1; -k -c];
B=[0; w];
disp = zeros(size(P));
x_zero = -inv(A)*B*P(1);
x_old = x_zero;
disp(1) = x_zero(1);
for pct = 1:length(P)-1
    x = inv(eye(size(A))-Delta_t*A)*(x_old + Delta_t*B*P(pct+1));
    disp(pct+1) = x(1);
    x_old = x;
end;
disp = disp - disp_shift;

%%% Extract the seven loops.

P_1 = P(141:281);
P_2 = P(406:526);
P_3 = P(631:731);
P_4 = P(816:896);
P_5 = P(961:1021);
P_6 = P(1066:1106);
P_7 = P(1131:1151);

disp_1 = disp(141:281);
disp_2 = disp(406:526);
disp_3 = disp(631:731);
disp_4 = disp(816:896);
disp_5 = disp(961:1021);
disp_6 = disp(1066:1106);
disp_7 = disp(1131:1151);

%%% Begin plots.

figure(1)
plot(Ec_gpts,nu_one,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([0 Ec_int 0 1.1*max(nu_one)]);
xlabel('V/m')
title('Coercive Field Density')

figure(2)
plot(Ei_gpts,nu_two,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-Ei_int Ei_int 0 1.1*max(nu_two)]);
xlabel('V/m')
title('Interactive Field Density')

figure(3)
plot(field_1,P_1*1e4,'k')
hold on
plot(field_2,P_2*1e4,'b')
plot(field_3,P_3*1e4,'r')
plot(field_4,P_4*1e4,'g')
plot(field_5,P_5*1e4,'c')
plot(field_6,P_6*1e4,'m')
plot(field_7,P_7*1e4,'y')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-1000 8000 0 1])
xlabel('Electric Field (V/m)');
ylabel('Polarization (C/m^2)');
hold off

figure(4)
plot(field_1,disp_1,'k')
hold on
plot(field_2,disp_2,'b')
plot(field_3,disp_3,'r')
plot(field_4,disp_4,'g')
plot(field_5,disp_5,'c')
plot(field_6,disp_6,'m')
plot(field_7,disp_7,'y')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-1000 8000 -5e-5 5e-5])
xlabel('Electric Field (V/m)');
ylabel('Displacement (m)');
hold off


