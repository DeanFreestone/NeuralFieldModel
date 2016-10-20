%% Numerical check. Implement equation (12)
% implementationf of Equation (12)

clc
clear
close all

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 201;

x = linspace(SpaceMin, SpaceMax, NPoints);

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x);

%% model parameters
% ~~~~~~~~~~~~~~~



tau = 0.01; % synaptic time constant

VtPlus1 = []; % field at T+1

% V at time t

%  v_t = ones(NPoints, NPoints);

rng(0,'twister');

lower = -6;
upper = 6;
Vt = (upper-lower).*rand(NPoints, NPoints) + lower;

% v_t = rand(NPoints, NPoints); % initialise a random field v at time point T



Ts = 0.0001; % time step

nx = 16; % number of Gaussian basis functions

theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

nTheta = 3; % number of connectivity kernel basis functions

% ~~~~~~~~~~~~~~~
% parameters for firing rate function
slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

ks = 1- Ts*(1/tau); % time constant parameter

errorPart = zeros(NPoints, NPoints); % set error part to zero for now
%% integral part


% initialise integral part
integralPart = zeros(NPoints, NPoints);


% firing rate function
firingRate_Vt = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - Vt)));


% integral. convolution or integral
for m = 1 : NPoints
    for n = 1 : NPoints
        r = [X(m, n), Y(m, n)]; % location r vector
        
        % connectivity kernel, a sum of Gaussian basis functions
        for p = 1 : nTheta
            gaussians(:,:, p) = Define2DGaussian_AnisotropicKernel(r(1), r(2), [vector_Sigma_Psi(p) 0; 0 vector_Sigma_Psi(p)], NPoints, SpaceMin, SpaceMax) * theta(p);
        end
        w = squeeze(sum(gaussians, 3)); % connectivity kernel
        
        % define connectivity kernel at location r
        integralPart = integralPart + w.*firingRate_Vt;
    end
end

%% v(t+1)

VtPlus1 = ks * Vt + Ts * integralPart + errorPart; % calculate v(t+1)
