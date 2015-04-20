function [outfield, varargout] = ferroic_hyst(infield, density, eta, varargin)
%
% function [outfield, varargout] = ferroic_hyst(infield, density, eta, varargin)       
%
% Computes the polarization or magnetization for a given electric or
% magnetic field.  The models are meant to accompany "Efficient Implementation 
% of Homogenized Free Energy Model with Emphasis on Thermal Relaxation"
% by Tom Braun and Ralph Smith, presented at the 2005 Smart Structures
% and Materials Conference.  This just implements the final optimized 
% models in each case.
%
% This core processing is written as a c language mex file, and the code
% must be compiled before this may be run.  To compile the file, type
% "mex ferroic_core.c" at the MATLAB prompt. 
%
% Inputs:
%     infield -- the input electric or magnetic field points
%     density -- structure with the following information
%         c_mode  -- Quadrature mode for coercive distribution integration.  This must be one of:
%             'trap' -- trapezoid rule
%             'simp' -- Simpson's rule
%             'simp38' -- Simpson's 3/8 rule
%             'gauss2' -- 2 point Gaussian quadrature
%             'gauss3' -- 3 point Gaussian quadrature
%             'gauss4' -- 4 point Gaussian quadrature
%         c_intervals -- number of intervals of the method c_mode to use for coercive field values
%         e_mode -- Quadrature mode for effective distribution integration.  For the relaxation
%             model this must be either 'trap' or 'simp38'
%         e_intervals -- number of intervals of the method e_mode to use for effective field values
%             this must give an even number of quadrature points when combined with c_mode
%         type -- 1 for normal/normal, 2 for lognormal (coercive)/normal(effective)
%                 3 for general densities.
%         If the type is 1, density should have the following other variables in the struct
%             C (a combination of c1 and c2, as given in the literature)
%             b_one
%             b_two
%         If the type is 2, it should have
%             C
%             c
%             Ec_bar
%             b_two
%         If the type is 3, it should have
%             Ec_int
%             Ee_int
%             nu_c -- a vector giving the density for Ec.  
%             nu_e -- a vector giving the density for Ee. 
%     eta -- a parameter to estimate -- permativity
%     beta -- relaxation model only.  A parameter to estimate (sqrt(2) k V / T)
%     tau -- relaxation model only.  A parameter to estimate
%     Delta_t -- relaxation model only.  Time between two samples of infield
%     resolution -- relaxation model only. An integer affecting the accuracy of the results
%         larger values will give more accurate results at the cost of more memory (and therefore often more time)
%     xpos -- optional: initial state of the material.  If not given, the material is assumed
%         start depoled.  Note that this can be a saved state from a previous run with the same parameters, 
%         thereby allowing the model to pick up where it left off.
% 
% LICENSE: This work is distributed under BSD license.  A complete copy
% of the license can be found at the top of model_core.c.  By running this
% code, you agreed to abide by the terms of the license.  
% 
% NOTE: This code is derived from a MATLAB homogenized free energy model 
% written by Andrew Hatch, North Carolina State University, Department 
% of Mathematics.  It is used with permission, and much thanks.  
% 
% Written by:
% Tom Braun
% North Carolina State University
% Department of Mathematics/Center for Research in Scientific Computation
% Raleigh NC  27695
% tbraun@pobox.com
% January 2005
%

if (nargin ~= 3 & nargin ~= 4 & nargin ~= 7 & nargin ~= 8)
    error('Wrong number of input arguments');
end
if (nargout > 2)
    error('More outputs were requested than the 2 that are available');
end

if density.type == 1
    c_one = sqrt(density.C);
    c_two = sqrt(density.C);
    Ec_int = sqrt(density.b_one*7);
    Ee_int = sqrt(density.b_two*7);
elseif density.type == 2
    c_one = sqrt(density.C);
    c_two = sqrt(density.C);
    Ec_int = density.Ec_bar*200^density.c;
    Ee_int = sqrt(density.b_two*7);
elseif density.type == 3
    Ec_int = density.Ec_int;
    Ee_int = density.Ee_int;
end

% Set up quadrature points and distributions for coercive field integration.
hc = Ec_int/density.c_intervals;
if (strcmp(density.c_mode, 'trap'))
    Nc = density.c_intervals+1;
    w_c = hc*ones(Nc,1);
    w_c(1) = hc/2;
    w_c(Nc) = hc/2;
    Ec_pts = linspace(0, Ec_int, Nc)';
elseif (strcmp(density.c_mode, 'simp'))
    Nc = 2*density.c_intervals + 1;
    w_c = hc*ones(Nc, 1)*2/3;
    w_c(1) = hc/6;
    w_c(Nc) = hc/6;
    v(3:2:Nc-2) = hc/3;
    Ec_pts = linspace(0, Ec_int, Nc)';
elseif (strcmp(density.c_mode, 'simp38'))
    Nc = 3*density.c_intervals + 1;
    w_c = hc*ones(Nc, 1)*9/24;
    w_c(1) = 3*hc/24;
    w_c(Nc) = 3*hc/24;
    w_c(4:3:Nc) = 3*hc/12;
    Ec_pts = linspace(0, Ec_int, Nc)';
elseif (strcmp(density.c_mode, 'gauss2'))
    Nc = 2*density.c_intervals;
    w_c = hc*ones(Nc, 1)/2;
    Ec_pts = zeros(Nc, 1);
    Ec_pts(1) = (1 - sqrt(1/3))*hc/2;
    Ec_pts(2) = (1 + sqrt(1/3))*hc/2;
    for index = 2:density.c_intervals;
        Ec_pts(2*index-1:2*index) = (index-1)*hc + Ec_pts(1:2);
    end
elseif (strcmp(density.c_mode, 'gauss3'))
    Nc = 3*density.c_intervals;
    w_c = hc*ones(Nc, 1)*0.27777777777777;
    w_c(2:3:Nc-1) = hc*0.44444444444444;
    Ec_pts = zeros(Nc, 1);
    Ec_pts(1) = (1/2 - 0.387298335)*hc;
    Ec_pts(2) = hc/2;
    Ec_pts(3) = (1/2 + 0.387298335)*hc;
    for index = 2:density.c_intervals;
        Ec_pts(3*index-2:3*index) = (index-1)*hc + Ec_pts(1:3);
    end
elseif (strcmp(density.c_mode, 'gauss4'))
    Nc = 4*density.c_intervals;
    w_c = hc*ones(Nc, 1)*49/(12*(18 + sqrt(30)));
    w_c(2:4:Nc) = hc*49/(12*(18 - sqrt(30)));
    w_c(3:4:Nc) = hc*49/(12*(18 - sqrt(30)));
    Ec_pts = zeros(Nc, 1);
    Ec_pts(1) = (1/2 - sqrt(15+2*sqrt(30))/(2*sqrt(35)))*hc;
    Ec_pts(2) = (1/2 - sqrt(15-2*sqrt(30))/(2*sqrt(35)))*hc;
    Ec_pts(3) = (1/2 + sqrt(15-2*sqrt(30))/(2*sqrt(35)))*hc;
    Ec_pts(4) = (1/2 + sqrt(15+2*sqrt(30))/(2*sqrt(35)))*hc;
    for index = 2:density.c_intervals;
        Ec_pts(4*index-3:4*index) = (index-1)*hc + Ec_pts(1:4);
    end
else
    error('Undefined mode for coercive field integration');
end
if density.type == 1
    nu_c = c_one*exp(-((Ec_pts).^2)/density.b_one);
elseif density.type == 2
    nu_c = c_one*exp(-(log(Ec_pts/density.Ec_bar)/(2*density.c)).^2);
elseif density.type == 3
    nu_c = density.nu_c;
end
w_c = w_c.*nu_c;

% Set up quadrature points and distributions for effective field distributions.
he = 2*Ee_int/density.e_intervals; % factor of two b/c this is neg-inf to pos-inf integrand
if (strcmp(density.e_mode, 'trap'))
    if nargin > 5 & mod(density.e_intervals, 2)~=1
        error('Algorithm 5 requires an odd number of density.e_intervals');
    end
    Ne = density.e_intervals+1;
    w_e = he*ones(Ne,1);
    w_e(1) = he/2;
    w_e(Ne) = he/2;
    Ee_pts = linspace(-Ee_int, Ee_int, Ne)';
elseif (strcmp(density.e_mode, 'simp'))
    if nargin > 5
        error('Algorithm 5 requires either trap or simp38 for density.e_mode');
    end
    Ne = 2*density.e_intervals + 1;
    w_e = he*ones(Ne, 1)*2/3;
    w_e(1) = he/6;
    w_e(Ne) = he/6;
    w_e(3:2:Ne-2) = he/3;
    Ee_pts = linspace(-Ee_int, Ee_int, Ne)';
elseif (strcmp(density.e_mode, 'simp38'))
    if nargin > 5 & mod(density.e_intervals, 2)~=1
        error('Algorithm 5 requires an odd number of density.e_intervals');
    end
    Ne = 3*density.e_intervals + 1;
    w_e = he*ones(Ne, 1)*9/24;
    w_e(1) = 3*he/24;
    w_e(Ne) = 3*he/24;
    w_e(4:3:Ne-3) = 3*he/12;
    Ee_pts = linspace(-Ee_int, Ee_int, Ne)';
elseif (strcmp(density.e_mode, 'gauss2'))
    if nargin > 5
        error('Algorithm 5 requires either trap or simp38 for density.e_mode');
    end
    Ne = 2*density.e_intervals;
    w_c = he*ones(Ne, 1)/2;
    Ee_pts = zeros(Ne, 1);
    Ee_pts(1) = (1 - sqrt(1/3))*he/2 - Ee_int;
    Ee_pts(2) = (1 + sqrt(1/3))*he/2 - Ee_int;
    for index = 2:density.e_intervals;
        Ee_pts(2*index-1:2*index) = (index-1)*he + Ee_pts(1:2);
    end
elseif (strcmp(density.e_mode, 'gauss3'))
    if nargin > 5
        error('Algorithm 5 requires either trap or simp38 for density.e_mode');
    end
    Ne = 3*density.e_intervals;
    w_e = hq*ones(Ne, 1)*0.27777777777777;
    w_e(2:3:Ne-1) = he*0.44444444444444;
    Ee_pts = zeros(Ne, 1);
    Ee_pts(1) = (1/2 - 0.387298335)*he - Ee_int;
    Ee_pts(2) = he/2 - Ee_int;
    Ee_pts(3) = (1/2 + 0.387298335)*he - Ee_int;
    for index = 2:density.e_intervals;
        Ee_pts(3*index-2:3*index) = (index-1)*he + Ee_pts(1:3);
    end
elseif (strcmp(density.e_mode, 'gauss4'))
    if nargin > 5
        error('Algorithm 5 requires either trap or simp38 for density.e_mode');
    end
    Ne = 4*density.e_intervals;
    w_e = he*ones(Ne, 1)*49/(12*(18 + sqrt(30)));
    w_e(2:4:Ne) = he*49/(12*(18 - sqrt(30)));
    w_e(3:4:Ne) = he*49/(12*(18 - sqrt(30)));
    Ee_pts = zeros(Ne, 1);
    Ee_pts(1) = (1/2 - sqrt(15+2*sqrt(30))/(2*sqrt(35)))*he - Ee_int;
    Ee_pts(2) = (1/2 - sqrt(15-2*sqrt(30))/(2*sqrt(35)))*he - Ee_int;
    Ee_pts(3) = (1/2 + sqrt(15-2*sqrt(30))/(2*sqrt(35)))*he - Ee_int;
    Ee_pts(4) = (1/2 + sqrt(15+2*sqrt(30))/(2*sqrt(35)))*he - Ee_int;
    for index = 2:density.e_intervals;
        Ee_pts(4*index-3:4*index) = (index-1)*he + Ee_pts(1:4);
    end
else
    error('Undefined mode for effective field integration');
end
if density.type == 1 | density.type == 2
    nu_e_pos = c_two*exp(-(Ee_pts(Ne/2+1:Ne).^2)/density.b_two);
elseif density.type == 3
    nu_e_pos = density.nu_e;
end
nu_e_neg = rot90(rot90(nu_e_pos));
nu_e = [nu_e_neg; nu_e_pos];
w_e = w_e.*nu_e;

if (nargin < 5)
    if (nargin == 3)
        xpos = ones(Ne, Nc);
        xpos(1:Ne / 2, :) = -1;
    else
        xpos = varargin{nargin - 3}
    end

    weightsum = sum(sum(w_c * w_e')) / eta;
    addit = sum(sum(w_c * (w_e .* Ee_pts)')) / eta;
    if (nargout == 2)
        [outfield, varargout{1}] = ferroic_core(infield, Ec_pts, w_c, Ee_pts, w_e, weightsum, addit, xpos);
    else
        outfield = ferroic_core(infield, Ec_pts, w_c, Ee_pts, w_e, weightsum, addit, xpos);
    end
else
    [beta, tau, Delta_t, resolution] = deal(varargin{1:4});
    if (nargin == 7)
        xpos = ones(Ne, Nc);
        xpos(1:Ne / 2, :) = 0;
    else
        xpos = varargin{nargin - 3}
    end

    weightsum = sum(sum(w_c * w_e'));
    addit = sum(sum(w_c * (w_e .* Ee_pts)')) / eta - weightsum;
    weightsum = weightsum / eta;
    epsilon = 1e-3 / beta;
    fudge = 1e-50;

    e_step = (Ee_pts(2) - Ee_pts(1)) / resolution;
    num_res = Ne * resolution;
    if (Ec_pts(Nc) > Ee_pts(Ne))
        increase = ceil((Ec_pts(Nc) - Ee_pts(Ne)) ./ e_step) + num_res;
    else
        increase = num_res - floor((Ee_pts(Ne) - Ec_pts(Nc)) ./ e_step);
    end
    step = e_step / (eta * beta);
    emu = [Ee_pts(1) ./ (eta * beta) - step * increase : step : Ee_pts(Ne) ./ (eta * beta) + step * increase];
    % To save memory, and because modern matlab versions will cache the value anyway,
    % we don't use temporary variables here as we did in C, but instead just include
    % some erfc functions calls twice.  This has been found to be quicker in Matlab
    % (but slower most everywhere else)
    p_negpos = (1/tau) .* ((sign(emu) + 1)/2 + (1 - sign(emu))/2 .* (1 - erfc(emu + epsilon) ./ (erfc(emu - epsilon) + fudge)));
    p_posneg = (1/tau) .* ((sign(emu) + 1)/2 .* (1 - erfc(-emu + epsilon) ./ (erfc(-emu - epsilon) + fudge)) + (1 - sign(emu))/2);
    P_neg = -(beta/sqrt(pi)) .* exp(-(emu + epsilon).^2) ./ (erfc(emu + epsilon) + fudge);
    P_pos = (beta/sqrt(pi)) .* exp(-(-emu + epsilon).^2) ./ (erfc(-emu + epsilon) + fudge);

    if (nargout == 2)
        [outfield, varargout{1}] = ferroic_core(infield, Ec_pts, w_c, e_step, w_e, weightsum, addit, xpos, p_negpos, p_posneg, P_neg, P_pos, ...
            Delta_t, increase, resolution);
    else
        outfield = ferroic_core(infield, Ec_pts, w_c, e_step, w_e, weightsum, addit, xpos, p_negpos, p_posneg, P_neg, P_pos, ...
            Delta_t, increase, resolution);
    end
end

