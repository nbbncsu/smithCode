%
%                      Program polarization_model_demo.m
%


clear all
%close all
format short e

%%% Set options.

freq_dep_yes = 0; % must be 0 or 1; 1 for model with thermal relaxation
density_type = 2; % 2 for lognormal/normal, 3 for general densities

%%% Load the model parameters.

if freq_dep_yes == 0
    if density_type == 2
        load params_pzt_log_demo.dat
        params = params_pzt_log_demo;
    elseif density_type == 3
        load params_pzt_gen_demo.dat
        params = params_pzt_gen_demo;
    end;
elseif freq_dep_yes == 1
    if density_type == 2
        load params_pzt_log_freq_demo.dat
        params = params_pzt_log_freq_demo;
    elseif density_type == 3
        load params_pzt_gen_freq_demo.dat
        params = params_pzt_gen_freq_demo;
    end;    
end;

%%% Set the value of Delta_t (which determines the relaxation).

Delta_t = 5e-2;
freq = 1/(100*Delta_t);

%%% Setup the input field.

E_1 = linspace(0,.5,26)';
E_2 = linspace(.5,-.5,51)';
E_3 = linspace(-.5,2,126)';
E_4 = linspace(2,-2,201)';
E_5 = linspace(-2,.4,121)';
E_6 = linspace(.4,-2,121)';
E_7 = linspace(-2,.6,131)';
E_8 = linspace(.6,-.3,46)';
E_9 = linspace(-.3,2,116)';
E_10 = linspace(2,-.7,136)';
E_11 = linspace(-.7,1,86)';
E_12 = linspace(1,-.5,76)';
E_13 = linspace(-.5,.5,51)';
E_14 = linspace(.5,-1.2,86)';
E_15 = linspace(-1.2,-.6,31)';
E_16 = linspace(-.6,-2,71)';

%%% The scaling factor at the end of the next command
%%% converts from units of MV/m to units of V/m.

E = [E_1(1:end-1); E_2(1:end-1); E_3(1:end-1); E_4(1:end-1);...
     E_5(1:end-1); E_6(1:end-1); E_7(1:end-1); E_8(1:end-1);...
     E_9(1:end-1); E_10(1:end-1); E_11(1:end-1); E_12(1:end-1);...
     E_13(1:end-1); E_14(1:end-1); E_15(1:end-1); E_16(1:end)]*1e6;

clear E_1 E_2 E_3 E_4 E_5 E_6 E_7 E_8 E_9 E_10 E_11 E_12 E_13 E_14 E_15 E_16
 
%%% Calculate the polarization given the input field.
 
[P, Ec_int, Ei_int, Ec_gpts, Ei_gpts, nu_one, nu_two] = ...
    polar_func(params, E, Delta_t, freq_dep_yes, density_type);

%%% Begin plots.

figure(1)
plot(Ec_gpts/1e6,nu_one,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([0 Ec_int/1e6 0 1.1*max(nu_one)]);
xlabel('MV/m')
title('Coercive Field Density')

figure(2)
plot(Ei_gpts/1e6,nu_two,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-Ei_int/1e6 Ei_int/1e6 0 1.1*max(nu_two)]);
xlabel('MV/m')
title('Interactive Field Density')

figure(3)
plot(E/1e6,P,'b')
hold on
plot([-3 3],[0 0],'k',[0 0],[-0.4 0.4],'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-3 3 -0.4 0.4])
xlabel('Electric Field (MV/m)');
ylabel('Polarization (C/m^2)');
hold off
