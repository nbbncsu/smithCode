function disp = rod_model(rodparams, field_in, Delta_t, model, modparams);

%
% function disp = rod_model(rodparams, field_in, Delta_t, model, modparams);
%
% Calculates the displacement of a ferroic rod using your choice of linear
% or nonlinear (homogenized energy model) magnetic field/magnetization (electic
% field/polarization) relations and either a lumped ordinary differential
% equation model that is simpler but only gives displacment at the end of 
% the rod, or a distributed partial differential equation model giving
% displacement along the whole rod.  Both the lumped and distributive model
% make use of trapezoid rule to numerically 'solve' the differential
% equation.
%
% Input arguments:
%
%   rodparams -- this is a structure containing all of the material dependent
%     parameters of the rod for the desired model.  
%     The elements in the structure are
%       youngmod  -- The young's modulus, in units of N/m^2
%       Cd        -- Kelvin-Voigt damping parameter for the rod, in Ns/m^2
%       density   -- density of the rod, in kg/m^3
%       rodlength -- length of the rod in meters
%       area      -- cross sectional area of the rod, in m^2
%       a         -- constant of proportionality relating magnetization to stress (in N)
%       ml        -- end mass attached to the free end of the rod, in kg
%       cl        -- end damping parameter of the rod, in Ns
%       kl        -- end stiffness of the rod, in N
% 
%   field_in -- the input electric or magnetic field, in A/m 
%     This should be a vector giving one input value for every period of
%     time.  Time periods are assumed to be equally spaced.  Note that
%     the rod is assumed to start at 0 displacement and depoled (in the 
%     case of the nonlinear homogenized energy model).  The input field
%     should take this into account.
%
%   Delta_t -- the time between two successive inputs (field_in).  This is
%     also the time between successive outputs
%
%   model -- the integer 1, 2, 3, or 4, where
%     1 is for the linear H-M (E-P), lumped ode model
%     2 is for the linear H-M (E-P), distributed pde model
%     3 is for the nonlinear homogenized energy H-M (E-P), lumped ode model
%     4 is for the nonlinear homogenized energy H-M (E-P), distributed pde model
%
%   modparams -- structure containing model parameters for the model specified 
%     in the previous parameter.  These will be different for different models.  
%     If running a linear H-M model, the following parameter should be in
%     the struct:
%       chi -- susceptibility (i.e. M = chi * H), which is unitless
%     If running the nonlinear homogenized energy H-M model, the following 
%     parameters should be in the struct:
%       density -- a struct containt the coercive and interactive field
%         distributions
%       eta
%       beta
%       tau
%     See the help for ferroic_hyst.m for details on these parameters.
%     For the lumpded ode model, no additional parameters are needed.  For
%     the distributed pde model, however, one more must also be specified:
%       N -- number of finite elements to use in the discretization.  This
%         should be an integer >= 2.
%
% Output
%   
%   disp -- the displacement of the rod.  If running the lumped ode model,
%     this will be vector with the same number of elements as field_in.  If
%     running the distributed pde model, this will by a matrix where with
%     the number of columns equal to the number of inputs and the number of 
%     rows is the number of nodes (modparams.N).  These are equally spaced from
%     rodlength/N to rodlength, with the displacement at 0 assumed to be 
%     zero for all time.
%
% Example 1: Linear E-P Model
%   rodparams.youngmod = 110e9;
%   rodparams.Cd = 3e6;
%   rodparams.density = 9250;
%   rodparams.rodlength = 0.115;
%   rodparams.area = 1.27e-4
%   rodparams.a = 2.75e-2;
%   rodparams.ml = 0.5;
%   rodparams.cl = 1e3;
%   rodparams.kl = 2e6;
%   Delta_t = 1e-3;
%   H = 5e4 * sin(2*pi*Delta_t*50*[0:5000]);
%   modparams.chi = 11.3;
%   modparams.N = 16;
%
% Example 2: Hysteretic E-P Model
%   rodparams.youngmod = 1.1e11;
%   rodparams.ml = 0.5;
%   rodparams.kl = 2e6;
%   rodparams.cl = 1000;
%   rodparams.Cd = 3e6;
%   rodparams.rodlength = 0.115;
%   rodparams.density = 9250;
%   rodparams.area = 1.2668e-4;
%   rodparams.a = 0.0275;
%   modparams.chi = 11.3;
%   modparams.N = 128;
%   modparams.eta = 6.9861e4;
%   modparams.beta = 2.3292e-5;
%   modparams.tau = 0.0690;
%   modparams.density.type = 1;
%   modparams.density.c_intervals = 21;
%   modparams.density.c_mode = 'simp38';
%   modparams.density.e_intervals = 21;
%   modparams.density.e_mode = 'simp38';
%   modparams.density.b_one = 6.0944e6;
%   modparams.density.b_two = 3.4814e6;
%   modparams.density.C = 0.0546;
%
% 
% Written by:
%
%   Tom Braun
%   North Carolina State University
%   Department of Mathematics/Center for Research in Scientific Computing
%   Raleigh, NC 27606
%
%   with the help of Ralph C. Smith
%   North Carolina State University
%   Department of Mathematics/Center for Research in Scientific Computing
%   Raleigh, NC 27695
%
% Copyright 2005 North Carolina State University.  Released under the terms
% of the BSD license, which can be found online at
% http://www.opensource.org/licenses/bsd-license.php, and which should have
% been distributed with this file.
% 

%%% Error checking  for rodparams (modparams will be checked later, once we know what model we're running) %%%
if (nargin ~= 5) error('wrong number of input arguments given to rod_model'); end
if (0==isfield(rodparams, 'youngmod')) error('rodparams must contain youngmod, the Young''s modulus'); end
if (0==isfield(rodparams, 'Cd')) error('rodparams must contain Cd, the Kelvin-Voigt damping parameter'); end
if (0==isfield(rodparams, 'density')) error('rodparams must contain density, the density of the rod'); end
if (0==isfield(rodparams, 'rodlength')) error('rodparams must contain rodlength, the length of the rod'); end
if (0==isfield(rodparams, 'area')) error('rodparams must contain area, the cross-sectional area of the rod'); end
if (0==isfield(rodparams, 'a')) error('rodparams must contain a, a constant relating magnetization/polarization to stress'); end
if (0==isfield(rodparams, 'ml')) error('rodparams must contain ml, the mass attached to the end of the rod'); end
if (0==isfield(rodparams, 'cl')) error('rodparams must contain cl, the damping on the end of the rod'); end
if (0==isfield(rodparams, 'kl')) error('rodparams must contain kl, the stiffness effecting the end of the rod'); end

%%% Convert magnetic or electric field to magnetization or polarization

if (model == 1 | model == 2)
    if (0==isfield(modparams, 'chi')) error('modparams must contain chi, the susceptiblity'); end
    magnetization = modparams.chi * field_in;
elseif (model == 3 | model == 4)
    magnetization = ferroic_hyst(field_in, modparams.density, modparams.eta, modparams.beta, modparams.tau, Delta_t, 10);
else
    error('Unrecognized model number given to rod model.  model must be 1, 2, 3, or 4');
end
if (model == 2 | model == 4)
    if (0==isfield(modparams, 'N')) error('modparams must contain N, the number of elements to discretize the rod into'); end
end

%%% Rod Model itself.  Convert magnetization into displacements %%%

% The lumped ode model is a scalar spring model.  The distributed model has
% the form of a PDE which is subsequently discretized in space using a linear
% finite element basis. 

if (model == 1 | model == 3) % lumped ode
    N = 1;
    M = rodparams.ml + rodparams.density * rodparams.area * rodparams.rodlength;
    C = rodparams.cl + rodparams.Cd * rodparams.area / rodparams.rodlength;
    K = rodparams.kl + rodparams.youngmod * rodparams.area / rodparams.rodlength;
    alpha = rodparams.a;
else % distributed pde
    % The integrals here are super simple, and have already been solved
    % when forming the matrices below, rather than performing Gaussian
    % quadrature each time the function is called.
    N = modparams.N;
    h = rodparams.rodlength ./ N;
    M = diag(ones(N,1) * 2 * h * rodparams.density * rodparams.area / 3, 0);
    M(N, N) = h * rodparams.density * rodparams.area / 3 + rodparams.ml;
    M = M + diag(ones(N-1,1) * h * rodparams.density * rodparams.area / 6, 1);
    M = M + diag(ones(N-1,1) * h * rodparams.density * rodparams.area / 6, -1);
    C = diag(ones(N,1) * 2 * rodparams.Cd * rodparams.area / h, 0);
    C(N, N) = rodparams.Cd * rodparams.area / h + rodparams.cl;
    C = C - diag(ones(N-1,1) * rodparams.Cd * rodparams.area / h, 1);
    C = C - diag(ones(N-1,1) * rodparams.Cd * rodparams.area / h, -1);
    K = diag(ones(N,1) * 2 * rodparams.youngmod * rodparams.area / h, 0);
    K(N, N) = rodparams.youngmod * rodparams.area / h + rodparams.kl;
    K = K - diag(ones(N-1,1) * rodparams.youngmod * rodparams.area / h, 1);
    K = K - diag(ones(N-1,1) * rodparams.youngmod * rodparams.area / h, -1);
    alpha = zeros(N, 1);
    alpha(N) = rodparams.a;
end

sysA=[zeros(N,N), eye(N, N); -inv(M)*K, -inv(M)*C];
sysB=[zeros(N,1); inv(M)*alpha];
tmp = inv(eye(size(sysA)) - Delta_t/2 * sysA);
sysA = tmp*(eye(size(sysA)) + Delta_t/2 * sysA);
sysB = Delta_t/2 * tmp * sysB;
disp = zeros(N, length(magnetization));
uk = zeros(2*N, 1);
for (index = 2:length(magnetization))
    uk = sysA * uk + sysB*(magnetization(index-1) + magnetization(index));
    disp(:, index) = uk(1:N);
end

