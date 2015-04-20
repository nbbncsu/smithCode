%
%         Rod_model_demo.m
%
%
%  Documentation for this code can be found in the function rod_model_construct.m 
%

  clear all
  close all

%
% Set parameters
%

  rodparams.youngmod = 110e9;
  rodparams.Cd = 3e6;
  rodparams.density = 9250;
  rodparams.rodlength = 0.115;
  rodparams.area = 1.27e-4
  rodparams.a = 2.75e-2;
  rodparams.ml = 0.5;
  rodparams.cl = 1e3;
  rodparams.kl = 2e6;
  Delta_t = 1e-4;
  time_int = 0:Delta_t:5;
  H = 5e4 * sin(2*pi*time_int);
  modparams.chi = 11.3;
  modparams.N = 32;

%
% Run model and plot solution
%
  disp = rod_model_construct(rodparams, H, Delta_t, 2, modparams);
  [Nd,td] = size(disp);
  end_disp = disp(Nd,:);

  plot(time_int, end_disp);

%
% End of Rod_model_demo.m
%

