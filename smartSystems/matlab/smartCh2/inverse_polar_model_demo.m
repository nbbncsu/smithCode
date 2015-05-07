clear all
%close all
format short e

%%% Set options.

freq_dep_yes = 0; % must be 0 (no thermal activation)
density_type = 2; % 2 for lognormal/normal, 3 for general densities

%%% Load the model parameters

if freq_dep_yes == 0
    if density_type == 2
        load params_pzt_log_demo.dat
        params = params_pzt_log_demo;
    elseif density_type == 3
        load params_pzt_gen_demo.dat
        params = params_pzt_gen_demo;
    end;
end;

%%% Set the value of Delta_t (which determines the frequency).

Delta_t = 5e-2;
freq = 1/(100*Delta_t);

%%% Setup the desired polarization.

P_1 = linspace(0,.025,6)';
P_2 = linspace(.025,-.025,11)';
P_3 = linspace(-.025,.35,76)';
P_4 = linspace(.35,-.35,141)';
P_5 = linspace(-.35,-.25,21)';
P_6 = linspace(-.25,-.35,21)';
P_7 = linspace(-.35,-.15,41)';
P_8 = linspace(-.15,-.2,11)';
P_9 = linspace(-.2,.35,111)';
P_10 = linspace(.35,.1,51)';
P_11 = linspace(.1,.3,41)';
P_12 = linspace(.3,.2,21)';
P_13 = linspace(.2,.25,11)';
P_14 = linspace(.25,-.25,101)';
P_15 = linspace(-.25,-.2,11)';
P_16 = linspace(-.2,-.35,31)';

P = [P_1(1:end-1); P_2(1:end-1); P_3(1:end-1);...
     P_4(1:end-1); P_5(1:end-1); P_6(1:end-1);...
     P_7(1:end-1); P_8(1:end-1); P_9(1:end-1);...
     P_10(1:end-1); P_11(1:end-1); P_12(1:end-1);...
     P_13(1:end-1); P_14(1:end-1); P_15(1:end-1);...
     P_16(1:end)];

clear P_1 P_2 P_3 P_4 P_5 P_6 P_7 P_8 P_9 P_10 P_11 P_12 P_13 P_14 P_15 P_16

%%% Create the time vector.

time = linspace(0,(length(P)-1)*Delta_t,length(P));

%%% Calculate the electric field given the polarization.

[E, P_out, Ec_int, Ei_int, Ec_gpts, Ei_gpts, nu_one, nu_two] = ...
    inverse_polar_func(params, P, Delta_t, freq_dep_yes, density_type);

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
plot(P,E/1e6,'b')
hold on
plot([-0.4 0.4],[0 0],'k',[0 0],[-2.1 2.1],'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-0.4 0.4 -2.1 2.1])
xlabel('Polarization (C/m^2)');
ylabel('Electric Field (MV/m)');
hold off

figure(4)
plot(E/1e6,P_out,'b')
hold on
plot([-2.1 2.1],[0 0],'k',[0 0],[-0.4 0.4],'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([-2.1 2.1 -0.4 0.4])
xlabel('Electric Field (MV/m)');
ylabel('Polarization (C/m^2)');
hold off

figure(5)
plot(time,P,'g')
hold on
plot(time,P_out,'r--')
plot([0 time(end)],[0 0],'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([0 time(end) -0.4 0.4])
xlabel('Time (s)');
ylabel('Polarization (C/m^2)');
legend('P_{in}','P_{out}',3)
hold off

figure(6)
plot(time, P - P_out,'k')
h = gca;
get(h);
set(h,'FontSize',[16]);
axis([0 time(end) 1.1*min(P-P_out) 1.1*max(P-P_out)])
xlabel('Time (s)');
ylabel('Polarization (C/m^2)');
hold off
