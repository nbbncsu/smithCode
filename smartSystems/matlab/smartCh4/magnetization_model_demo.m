%
%                    Magnetization_model_demo.m
%

clear all
%close all
format short e

%%% Set options.

freq_dep_yes = 0; % must be 0 or 1, 1 for thermal activation
density_type = 1; % 1 for normal/normal density choice

%%% Load the model parameters.

if freq_dep_yes == 0
    if density_type == 1
        load params_ter_nor_demo.dat
        params = params_ter_nor_demo;
    end;
elseif freq_dep_yes == 1
    if density_type == 1
        load params_ter_nor_freq_demo.dat
        params = params_ter_nor_freq_demo;
    end;    
end;

%%% Set the value of Delta_t (which determines the thermal activation level).

Delta_t = 5e-2;
freq = 1/(100*Delta_t);

%%% Setup the input field.

H_1 = linspace(0,.5,26)';
H_2 = linspace(.5,-.5,51)';
H_3 = linspace(-.5,2,126)';
H_4 = linspace(2,-2,201)';
H_5 = linspace(-2,.4,121)';
H_6 = linspace(.4,-2,121)';
H_7 = linspace(-2,.6,131)';
H_8 = linspace(.6,-.3,46)';
H_9 = linspace(-.3,2,116)';
H_10 = linspace(2,-.7,136)';
H_11 = linspace(-.7,1,86)';
H_12 = linspace(1,-.5,76)';
H_13 = linspace(-.5,.5,51)';
H_14 = linspace(.5,-1.2,86)';
H_15 = linspace(-1.2,-.6,31)';
H_16 = linspace(-.6,-2,71)';

%%% The scaling factor at the end of the next command
%%% gives the magnetic field the correct order of
%%% magnitude where the units are A/m.

H = [H_1(1:end-1); H_2(1:end-1); H_3(1:end-1); H_4(1:end-1);...
     H_5(1:end-1); H_6(1:end-1); H_7(1:end-1); H_8(1:end-1);...
     H_9(1:end-1); H_10(1:end-1); H_11(1:end-1); H_12(1:end-1);...
     H_13(1:end-1); H_14(1:end-1); H_15(1:end-1); H_16(1:end)]*3e4;
 
clear H_1 H_2 H_3 H_4 H_5 H_6 H_7 H_8 H_9 H_10 H_11 H_12 H_13 H_14 H_15 H_16

%%% Calculate the magnetization given the magnetic field.
%%% The function mag_func is the same as polar_func used in the
%%% The function polar_func is called even though we are using
%%% polarization code since the polarization and magnetization
%%% models are identical from a programming perspective.
 
[M, Hc_int, Hi_int, Hc_gpts, Hi_gpts, nu_one, nu_two] = ...
    mag_func(params, H, Delta_t, freq_dep_yes, density_type);

%%% Begin plots.

figure(1)
plot(Hc_gpts,nu_one,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([0 Hc_int 0 1.1*max(nu_one)]);
xlabel('A/m')
title('Coercive Field Density')

figure(2)
plot(Hi_gpts,nu_two,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-Hi_int Hi_int 0 1.1*max(nu_two)]);
xlabel('A/m')
title('Interactive Field Density')

figure(3)
plot(H,M,'b')
hold on
plot([-6e4 6e4],[0 0],'k',[0 0],[-5e5 5e5],'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-6e4 6e4 -5e5 5e5])
xlabel('Magnetic Field (A/m)');
ylabel('Magnetization (A/m)');
hold off


